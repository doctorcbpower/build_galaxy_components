"""The GalIC + gas-disc workflow: report, snapshot checks, add-gas, command line, script and notebook."""
import ast
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from galaxy_components import cli, gas_disc, workflow
from galaxy_components.params import read_params
from galaxy_components.system import read_system

EX = Path(__file__).resolve().parents[1] / "examples" / "galic_gas_disc"
PARAMS = EX / "galic_m1.toml"


@pytest.fixture(scope="module")
def p():
    q = read_params(PARAMS)
    q.ndisc = 3000
    q.nrbin, q.nzbin = 60, 250
    return q


@pytest.fixture(scope="module")
def snapshot(p, tmp_path_factory):
    return workflow.standin_snapshot(p, tmp_path_factory.mktemp("snap") / "snap.hdf5", n_halo=4000, n_disc=4000, n_bulge=1000)


def test_report_rows_and_resolution_scaling(p):
    rep = workflow.model_report(p)
    assert [r["R_over_Rd"] for r in rep.rows] == list(workflow.RADII)
    assert all(r["vrot"] > 0 and r["H"] > 0 and r["h_over_H"] > 0 for r in rep.rows)
    big = workflow.model_report(p, ndisc=8 * p.ndisc)
    for a, b in zip(rep.rows, big.rows):
        assert b["h_over_H"] == pytest.approx(a["h_over_H"] / 2, rel=1e-6)
    assert "Vertical resolution" in rep.format() and "GalIC convention" in rep.format()
    # the model's halo is lighter than GalIC's by exactly the gas mass
    m = rep.model
    assert (m.m200 - m.mbulge - m.mstar - m.mbh) - m.mhalo == pytest.approx(m.mdisc - m.mstar)


def test_standin_matches_model(p, snapshot):
    snap = read_system(snapshot)
    finds = workflow.compare_snapshot(p, workflow.model_report(p).model, snap.species)
    assert len(finds) == 7 and all(f.ok for f in finds), [f.format() for f in finds if not f.ok]


def test_compare_flags_a_thin_model_against_a_thick_snapshot(p, snapshot):
    thin = read_params(PARAMS)
    thin.zdisc = 0.0
    model = workflow.model_report(thin).model
    finds = {f.name: f for f in workflow.compare_snapshot(thin, model, read_system(snapshot).species)}
    assert not finds["disc thickness z0 [kpc]"].ok


def test_compare_flags_reversed_spin(p, snapshot):
    snap = read_system(snapshot)
    snap.species[2]["vel"] = -snap.species[2]["vel"]
    finds = {f.name: f for f in workflow.compare_snapshot(p, workflow.model_report(p).model, snap.species)}
    assert not finds["disc spin direction (Lz > 0)"].ok


def test_add_gas_to_snapshot(p, snapshot, tmp_path):
    out = tmp_path / "ic.hdf5"
    with pytest.warns(UserWarning, match="under-resolved"):
        res = workflow.add_gas_to_snapshot(p, snapshot, out, seed=4, verbose=False)
    assert res["counts"] == {0: p.ndisc, 1: 4000, 2: 4000, 3: 1000}
    ic = read_system(out)
    rep = workflow.model_report(p)
    r0, rcut, m0, mkeep = gas_disc.disc_window(p, rep.disc)
    assert ic.mass(0) == pytest.approx(mkeep - m0, rel=1e-5)
    assert ic.mass(1) == pytest.approx(read_system(snapshot).mass(1), rel=1e-6)   # snapshot untouched
    assert np.all(np.isfinite(ic.species[0]["ug"]))
    # same seed -> same gas
    out2 = tmp_path / "ic2.hdf5"
    with pytest.warns(UserWarning):
        workflow.add_gas_to_snapshot(p, snapshot, out2, seed=4, verbose=False)
    assert np.array_equal(read_system(out2).species[0]["pos"], ic.species[0]["pos"])


def test_add_gas_refuses_gas_and_strict_mismatch(p, snapshot, tmp_path):
    gas = tmp_path / "has_gas.hdf5"
    with pytest.warns(UserWarning):
        workflow.add_gas_to_snapshot(p, snapshot, gas, verbose=False)
    with pytest.raises(ValueError, match="already contains gas"):
        workflow.add_gas_to_snapshot(p, gas, tmp_path / "x.hdf5", verbose=False)
    thin = read_params(PARAMS)
    thin.zdisc, thin.ndisc = 0.0, 500
    with pytest.raises(ValueError, match="does not match"):
        workflow.add_gas_to_snapshot(thin, snapshot, tmp_path / "y.hdf5", strict=True, verbose=False)


def test_cli_check_and_add_gas(snapshot, tmp_path, capsys):
    cli.main(["check", str(PARAMS), "--ndisc", "5000"])
    assert "ndisc = 5,000" in capsys.readouterr().out
    out = tmp_path / "cli_ic.hdf5"
    with pytest.warns(UserWarning):
        cli.main(["add-gas", str(PARAMS), "--snapshot", str(snapshot), "--out", str(out), "--ndisc", "800",
                  "--sampling", "poisson", "--strict"])
    assert read_system(out).counts()[0] == 800


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_build_script(snapshot, tmp_path):
    out = tmp_path / "script_ic.hdf5"
    env = {"GC": f"{sys.executable} -W ignore -m galaxy_components", "NDISC": "500", "STRICT": "1",
           "PYTHONPATH": str(Path(__file__).resolve().parents[1]), "PATH": "/usr/bin:/bin"}
    r = subprocess.run([str(EX / "build.sh"), str(snapshot), str(out)], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    assert read_system(out).counts()[0] == 500


def test_notebook_is_valid_and_its_code_parses():
    nb = json.loads((EX / "build_gas_disc_ics.ipynb").read_text())
    assert nb["nbformat"] == 4
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    assert len(code) >= 8
    for c in code:
        ast.parse("".join(c["source"]))          # no syntax errors, no magics
    assert all(o.get("output_type") != "error" for c in code for o in c.get("outputs", []))   # saved from a clean run


def test_notebook_executes(tmp_path, monkeypatch):
    nbformat = pytest.importorskip("nbformat")
    nbclient = pytest.importorskip("nbclient")
    pytest.importorskip("matplotlib")
    monkeypatch.setenv("PYTHONPATH", str(Path(__file__).resolve().parents[1]))   # the kernel imports the checkout
    work = tmp_path / "nb"
    shutil.copytree(EX, work)
    nb = nbformat.read(work / "build_gas_disc_ics.ipynb", as_version=4)
    for c in nb.cells:                           # keep it quick
        if c.cell_type == "code":
            c.source = c.source.replace("NDISC = 100_000", "NDISC = 4000")
    nbclient.NotebookClient(nb, timeout=600, resources={"metadata": {"path": str(work)}}).execute()
    assert (work / "ic_with_gas.hdf5").exists()
