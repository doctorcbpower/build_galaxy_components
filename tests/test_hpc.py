"""hpc/setup_run.py: run directories, SLURM scripts and the Arepo parameter file written from the ICs."""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

from galaxy_components import workflow
from galaxy_components.params import read_params
from galaxy_components.system import read_system

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "hpc" / "cases"
spec = importlib.util.spec_from_file_location("setup_run", ROOT / "hpc" / "setup_run.py")
setup_run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(setup_run)


@pytest.fixture(scope="module")
def shim(tmp_path_factory):
    """`galaxy-components` and `python3` on PATH that use this checkout's interpreter and sources."""
    d = tmp_path_factory.mktemp("bin")
    (d / "galaxy-components").write_text(f'#!/bin/bash\nexec env PYTHONPATH={ROOT} {sys.executable} -W ignore -m galaxy_components "$@"\n')
    (d / "python3").write_text(f'#!/bin/bash\nexec env PYTHONPATH={ROOT} {sys.executable} "$@"\n')
    for f in d.iterdir():
        f.chmod(0o755)
    return d


@pytest.fixture()
def cluster(tmp_path):
    f = tmp_path / "cluster.env"
    f.write_text('ACCOUNT=proj1   # a comment\nPARTITION=\nMODULES="gcc openmpi"\nAREPO_BIN=/bin/true\n'
                 'AREPO_NODES=1\nAREPO_NTASKS_PER_NODE=4\nAREPO_TIME=02:00:00\nENV_SETUP=source /x/bin/activate\n')
    return f


def make(case, tmp_path, cluster, scale=0.002, run=None):
    d = tmp_path / case
    setup_run.main(["--galaxy", str(CASES / case / "galaxy.toml"), "--run", str(run or CASES / case / "run.toml"),
                    "--cluster", str(cluster), "--run-dir", str(d), "--scale", str(scale)])
    return d


def run_ics(d, shim, nomodule=True):
    text = (d / "make_ics.sbatch").read_text()
    if nomodule:                                   # no `module` command outside a cluster
        text = "\n".join(l for l in text.splitlines() if not l.startswith(("module ", "source ")))
        (d / "make_ics_local.sh").write_text(text + "\n")
    r = subprocess.run(["bash", str(d / "make_ics_local.sh")], capture_output=True, text=True,
                       env=dict(os.environ, PATH=f"{shim}:{os.environ['PATH']}"))
    assert r.returncode == 0, r.stdout[-600:] + r.stderr[-800:]


def test_cluster_file_parsing(cluster):
    cl = setup_run.read_env(cluster)
    assert cl["ACCOUNT"] == "proj1" and cl["PARTITION"] == "" and cl["MODULES"] == "gcc openmpi"
    assert cl["AREPO_BIN"] == "/bin/true" and cl["LAUNCHER"] == "srun" and cl["ENV_SETUP"] == "source /x/bin/activate"


def test_run_file_is_strict(tmp_path):
    bad = tmp_path / "run.toml"
    bad.write_text("[run]\nt_end = 5\n")
    with pytest.raises(SystemExit, match="unknown field 't_end'"):
        setup_run.read_run(bad)
    bad.write_text("[runn]\nname = 'x'\n")
    with pytest.raises(SystemExit, match="unknown table"):
        setup_run.read_run(bad)


def test_empty_sbatch_values_drop_the_line():
    out = setup_run.render("#SBATCH --account=@A@\n#SBATCH --time=@T@\necho @T@\n", dict(A="", T="1:00"))
    assert out == "#SBATCH --time=1:00\necho 1:00\n"
    with pytest.raises(SystemExit, match="no value"):
        setup_run.render("echo @MISSING@\n", {})


@pytest.mark.parametrize("case", ["hot_halo", "hot_halo_disc"])
def test_cases_build_ics_and_arepo_parameters(case, tmp_path, cluster, shim):
    d = make(case, tmp_path, cluster)
    ics, arepo = (d / "make_ics.sbatch").read_text(), (d / "arepo.sbatch").read_text()
    assert "#SBATCH --account=proj1" in ics and "--partition" not in ics and "module load gcc" in ics
    assert "--background-grid" in ics and "--ntasks-per-node=4" in arepo and "--mem-per-cpu" not in arepo
    for script in ("make_ics.sbatch", "arepo.sbatch", "submit.sh"):
        assert subprocess.run(["bash", "-n", str(d / script)]).returncode == 0
    run_ics(d, shim)
    ic = read_system(d / "ICs" / "ic.hdf5")
    assert ic.counts()[0] > 0 and ic.counts()[1] > 0
    par = (d / "arepo_param.txt").read_text()
    assert "@" not in par.replace("%  ", "")        # every placeholder resolved
    hdr = setup_run.ic_header(d / "ICs" / "ic.hdf5")
    assert hdr["boxsize"] > 0 and f"BoxSize                       {hdr['boxsize']:.8g}" in par
    # the gas mass the cells are refined towards is the particle mass, not the mean (unequal-mass) cell mass
    import json
    sug = json.loads((d / "ICs" / "ic.hdf5.arepo.json").read_text())
    assert hdr["gas_mass"] == pytest.approx(sug["ReferenceGasPartMass"]) and f"{sug['ReferenceGasPartMass']:.8g}" in par
    assert hdr["gas_mass"] < ic.mass(0) / ic.counts()[0]
    assert "TimeMax                       5.11356" in par and "TimeLimitCPU                  6300" in par
    if case == "hot_halo_disc":                      # a cold disc: gas particles near the plane inside the halo gas
        assert ic.mass(0) > 0


def test_galic_snapshot_case(tmp_path, cluster, shim):
    p = read_params(CASES / "hot_halo_disc_galic" / "galaxy.toml")
    snap = workflow.standin_snapshot(p, tmp_path / "galic.hdf5", n_halo=3000, n_disc=2000, n_bulge=500)
    run = tmp_path / "run.toml"
    run.write_text((CASES / "hot_halo_disc_galic" / "run.toml").read_text().replace('snapshot = ""', f'snapshot = "{snap}"'))
    d = make("hot_halo_disc_galic", tmp_path, cluster, run=run)
    assert "merge merge.toml" in (d / "make_ics.sbatch").read_text() and "nbody_components = false" in (d / "merge.toml").read_text()
    run_ics(d, shim)
    ic = read_system(d / "ICs" / "ic.hdf5")
    assert [ic.counts()[t] for t in (1, 2, 3)] == [3000, 2000, 500]
    assert ic.counts()[0] > 0 and (d / "arepo_param.txt").exists()
    # the collisionless particles are the snapshot's, untouched in mass (the box shift moves positions only)
    assert ic.mass(1) == pytest.approx(read_system(snap).mass(1), rel=1e-6)


def test_missing_snapshot_is_an_error(tmp_path, cluster):
    with pytest.raises(SystemExit, match="no N-body particles"):
        make("hot_halo_disc_galic", tmp_path, cluster)            # no snapshot and no N-body components: no dark matter
    run = tmp_path / "run.toml"
    run.write_text((CASES / "hot_halo_disc_galic" / "run.toml").read_text().replace('snapshot = ""', 'snapshot = "/nope/galic.hdf5"'))
    with pytest.raises(SystemExit, match="does not exist"):
        make("hot_halo_disc_galic", tmp_path, cluster, run=run)


def test_arepo_job_refuses_without_binary_or_ics(tmp_path, cluster):
    cluster.write_text("AREPO_BIN=\n")
    d = make("hot_halo", tmp_path, cluster)
    r = subprocess.run(["bash", str(d / "arepo.sbatch")], capture_output=True, text=True)
    assert r.returncode == 1 and "AREPO_BIN" in r.stderr
    cluster.write_text(f"AREPO_BIN={sys.executable}\n")
    d = make("hot_halo", tmp_path / "again", cluster)
    r = subprocess.run(["bash", str(d / "arepo.sbatch")], capture_output=True, text=True)
    assert r.returncode == 1 and "make_ics.sbatch first" in r.stderr
