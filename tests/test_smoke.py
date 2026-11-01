"""Smoke tests: run every command end to end on tiny inputs and sanity-check what is written.

These do not check the physics (the other test files do); they catch a command that
crashes, writes nothing, or writes NaNs / inconsistent files.   pytest -m smoke
"""
import subprocess
import sys
import tomllib
from pathlib import Path

import h5py
import numpy as np
import pytest

from galaxy_components.cli import main

pytestmark = pytest.mark.smoke

ROOT = Path(__file__).parent.parent
DATA = Path(__file__).parent / "data"

TINY = """
&MODEL
nrbin=40
nzbin=200
/
"""


@pytest.fixture(scope="module")
def params(tmp_path_factory):
    """A small parameter file: coarse disc mesh, few particles."""
    d = tmp_path_factory.mktemp("smoke")
    text = (DATA / "galaxy_params.txt").read_text()
    for old, new in (("ndisc=20000", "ndisc=600"), ("nhalo=20000", "nhalo=800"),
                     ("npart_comp(1)=20000", "npart_comp(1)=1500")):
        assert old in text, old
        text = text.replace(old, new)
    path = d / "p.txt"
    path.write_text(text + TINY)
    return path


def check_ic(path, expect_types=None, boxed=False):
    """The file exists, is internally consistent, and contains no junk."""
    path = Path(path)
    assert path.exists() and path.stat().st_size > 0, f"{path} not written"
    with h5py.File(path) as f:
        h = f["Header"].attrs
        counts = np.asarray(h["NumPart_ThisFile"])
        masstab = np.asarray(h["MassTable"])
        ids = []
        for t in range(6):
            key = f"PartType{t}"
            if counts[t] == 0:
                assert key not in f
                continue
            g = f[key]
            n = counts[t]
            assert g["Coordinates"].shape == (n, 3) and g["Velocities"].shape == (n, 3)
            for name in ("Coordinates", "Velocities", "ParticleIDs"):
                assert np.all(np.isfinite(g[name][:])), f"{key}/{name} has non-finite values"
            if "Masses" in g:
                m = g["Masses"][:]
                assert m.shape == (n,) and np.all(np.isfinite(m)) and m.sum() > 0
                # zero-mass gas cells are the point of the AREPO background grid; nothing else may be massless
                assert np.all(m >= 0) if (t == 0 and boxed) else np.all(m > 0)
            else:
                assert masstab[t] > 0, f"{key}: no Masses block and no MassTable entry"
            if t == 0:
                assert "InternalEnergy" in g and np.all(g["InternalEnergy"][:] >= 0)
                assert np.all(np.isfinite(g["InternalEnergy"][:]))
            ids.append(g["ParticleIDs"][:])
            if boxed:
                x = g["Coordinates"][:]
                assert np.all((x > 0) & (x < h["BoxSize"])), f"{key} outside the box"
        ids = np.concatenate(ids)
        assert len(np.unique(ids)) == len(ids), "duplicate particle IDs"
        if boxed:
            assert h["BoxSize"] > 0
        if expect_types is not None:
            present = {t for t in range(6) if counts[t] > 0}
            assert present == set(expect_types), f"types {present}, expected {set(expect_types)}"
    return counts


# ------------------------------------------------------------------ single-component commands
def test_disc(params, tmp_path):
    main(["disc", str(params), "--outdir", str(tmp_path)])
    check_ic(tmp_path / "gas_disc.gdt.hdf5", {0})
    for name in ("vrot_components.txt", "disc_parameters.txt"):
        assert np.all(np.isfinite(np.loadtxt(tmp_path / name)[:, :3]))


def test_nbody(params, tmp_path):
    main(["nbody", str(params), "--outdir", str(tmp_path)])
    check_ic(tmp_path / "nfw_halo.gdt.hdf5", {1})
    assert (tmp_path / "distribution_function.txt").stat().st_size > 0


@pytest.mark.parametrize("flags", [[], ["--legacy", "--cold"], ["--halo-core", "0.1"], ["--zero-edge"]])
def test_halo(params, tmp_path, flags):
    main(["halo", str(params), "--outdir", str(tmp_path), *flags])
    check_ic(tmp_path / "gas_halo.gdt.hdf5", {0})


# ------------------------------------------------------------------ whole systems
@pytest.mark.parametrize("flags,boxed", [
    ([], False),
    (["--background-grid"], True),
    (["--halo-core", "0.1", "--zero-edge", "--background-grid"], True),
    (["--naive-merge"], False),
    (["--legacy-halo", "--cold-halo", "--no-nbody"], False),
])
def test_galaxy_variants(params, tmp_path, flags, boxed):
    main(["galaxy", str(params), "--outdir", str(tmp_path), "--name", "g.hdf5", *flags])
    types = {0} if "--no-nbody" in flags else {0, 1}
    counts = check_ic(tmp_path / "g.hdf5", types, boxed=boxed)
    assert counts[0] > 0


def test_system_alias(params, tmp_path):
    main(["system", str(params), "--outdir", str(tmp_path), "--name", "s.hdf5"])
    check_ic(tmp_path / "s.hdf5", {0, 1})


# ------------------------------------------------------------------ merging
def _merge_config(tmp_path, params, extra=""):
    (tmp_path / "gal.txt").write_text(Path(params).read_text())
    cfg = tmp_path / "m.toml"
    cfg.write_text(f"""
output = "merger.hdf5"
seed = 3
{extra}
[orbit]
systems = [0, 1]
r_init = 800.0
r_peri = 30.0
ecc = 1.0
[[system]]
name = "a"
params = "gal.txt"
euler = [0, 30, 0]
[system.overrides]
ndisc = 300
nhalo = 400
npart_comp = [800, 0, 0, 0, 0]
[[system]]
name = "b"
params = "gal.txt"
type_map = {{ 1 = 2 }}
[system.options]
halo_core = 0.1
[system.overrides]
m200 = 40.0
ndisc = 200
nhalo = 300
npart_comp = [500, 0, 0, 0, 0]
""")
    return cfg


@pytest.mark.parametrize("extra,boxed", [
    ("[background_grid]\nenabled = true", True),
    ('[gas_overlap]\nmode = "none"', False),
])
def test_merge(params, tmp_path, extra, boxed):
    cfg = _merge_config(tmp_path, params, extra)
    main(["merge", str(cfg), "--outdir", str(tmp_path / "out")])
    check_ic(tmp_path / "out" / "merger.hdf5", {0, 1, 2}, boxed=boxed)


def test_example_config_is_valid():
    """The shipped example parses and refers to files that exist."""
    cfg = tomllib.loads((ROOT / "examples" / "merger.toml").read_text())
    assert len(cfg["system"]) == 2 and cfg["orbit"]["systems"] == [0, 1]
    for s in cfg["system"]:
        assert (ROOT / "examples" / s["params"]).exists()


# ------------------------------------------------------------------ packaging / entry point
def test_module_entry_point(params, tmp_path):
    r = subprocess.run([sys.executable, "-m", "galaxy_components", "--help"],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0
    for cmd in ("disc", "nbody", "halo", "galaxy", "system", "merge"):
        assert cmd in r.stdout
    r = subprocess.run([sys.executable, "-m", "galaxy_components", "nbody", str(params),
                        "--outdir", str(tmp_path)], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-500:]
    check_ic(tmp_path / "nfw_halo.gdt.hdf5", {1})


def test_bad_input_fails_cleanly(tmp_path):
    with pytest.raises(FileNotFoundError):
        main(["nbody", str(tmp_path / "nope.txt")])


def test_check_ic_is_not_vacuous(tmp_path):
    """The checker really rejects broken files (NaNs, duplicate IDs, missing masses)."""
    from galaxy_components.io import write_hdf5
    good = dict(pos=np.ones((4, 3)), vel=np.zeros((4, 3)), mass=1.0)
    write_hdf5(tmp_path / "ok.hdf5", {1: good}, verbose=False)
    check_ic(tmp_path / "ok.hdf5", {1})
    bad = dict(good, pos=np.full((4, 3), np.nan))
    write_hdf5(tmp_path / "nan.hdf5", {1: bad}, verbose=False)
    with pytest.raises(AssertionError):
        check_ic(tmp_path / "nan.hdf5")
    write_hdf5(tmp_path / "nomass.hdf5", {1: dict(good, mass=0.0)}, verbose=False)
    with pytest.raises(AssertionError):
        check_ic(tmp_path / "nomass.hdf5")
    with pytest.raises(AssertionError):
        check_ic(tmp_path / "ok.hdf5", {0})                      # wrong particle types
    with pytest.raises(AssertionError):
        check_ic(tmp_path / "absent.hdf5")
    gas_no_u = {0: good}
    write_hdf5(tmp_path / "gas.hdf5", gas_no_u, verbose=False)   # gas without InternalEnergy
    with pytest.raises(AssertionError):
        check_ic(tmp_path / "gas.hdf5")


# ------------------------------------------------------------------ the standalone module examples
EXAMPLES = [("01_gas_disc.py", True, {0}), ("02_gas_halo.py", False, None),
            ("03_nbody_halo.py", False, None), ("04_background_grid.py", True, {0}),
            ("05_stitch_systems.py", True, {0, 1, 2}), ("06_custom_potential.py", False, None),
            ("07_anisotropic_halo.py", False, None)]


@pytest.mark.parametrize("script,writes,types", EXAMPLES)
def test_examples_run(script, writes, types, tmp_path):
    """Every script in examples/ runs from the repo root and (where it writes) writes a valid IC."""
    import os
    cmd = [sys.executable, str(ROOT / "examples" / script)]
    out = tmp_path / "out.hdf5"
    if writes:
        cmd += ["--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT,
                       env={**os.environ, "PYTHONPATH": str(ROOT)})
    assert r.returncode == 0, r.stderr[-800:]
    assert r.stdout.strip(), "example printed nothing"
    assert "nan" not in r.stdout.lower()
    if writes:
        check_ic(out, types, boxed=script in ("04_background_grid.py", "05_stitch_systems.py"))


# ------------------------------------------------------------------ packaging metadata
def test_pyproject_is_consistent():
    """Version, dependencies and the console script in pyproject.toml match the code."""
    import importlib

    import galaxy_components
    meta = tomllib.loads((ROOT / "pyproject.toml").read_text())
    proj = meta["project"]
    assert proj["name"] == "galaxy-components" and "version" in proj["dynamic"]
    assert meta["tool"]["setuptools"]["dynamic"]["version"]["attr"] == "galaxy_components.__version__"
    assert isinstance(galaxy_components.__version__, str) and galaxy_components.__version__.count(".") == 2
    mod, func = proj["scripts"]["galaxy-components"].split(":")
    assert callable(getattr(importlib.import_module(mod), func))
    name = lambda d: d.split(">")[0].split("=")[0]
    deps = {name(d) for d in proj["dependencies"]}
    optional = {name(d) for ds in proj["optional-dependencies"].values() for d in ds}
    assert {"numpy", "scipy", "h5py"} <= deps and "f90nml" in optional and "f90nml" not in deps
    # every third-party module the package imports is declared
    imported = set()
    for f in (ROOT / "galaxy_components").glob("*.py"):
        for line in f.read_text().splitlines():
            line = line.strip()
            if line.startswith(("import ", "from ")) and not line.startswith("from ."):
                imported.add(line.split()[1].split(".")[0])
    third_party = imported - set(sys.stdlib_module_names) - {"galaxy_components"}
    assert third_party <= deps | optional, f"undeclared dependencies: {third_party - deps - optional}"
    # optional packages must be imported lazily (inside a function), never at module level
    for f in (ROOT / "galaxy_components").glob("*.py"):
        for line in f.read_text().splitlines():
            assert not any(line.startswith((f"import {o}", f"from {o}")) for o in optional), f"{f.name}: {line}"


# ------------------------------------------------------------------ TOML parameter files
def test_convert_then_run_from_toml(params, tmp_path):
    """convert (namelist -> TOML) writes a file every command can then use."""
    toml = tmp_path / "p.toml"
    main(["convert", str(params), "-o", str(toml)])
    assert toml.read_text().startswith("# Converted from")
    main(["nbody", str(toml), "--outdir", str(tmp_path)])
    check_ic(tmp_path / "nfw_halo.gdt.hdf5", {1})
    main(["galaxy", str(toml), "--outdir", str(tmp_path), "--name", "t.hdf5"])
    check_ic(tmp_path / "t.hdf5", {0, 1})


def test_convert_to_stdout(params, capsys):
    main(["convert", str(params)])
    out = capsys.readouterr().out
    assert "[halo]" in out and "m200" in out


def test_example_toml_runs_a_system(tmp_path):
    """The shipped TOML example is a working parameter file (shrunk so it runs quickly)."""
    from galaxy_components.system import build_system
    s = build_system(ROOT / "examples" / "galaxy_params.toml", "x",
                     overrides=dict(ndisc=300, nhalo=500, npart_comp=[800, 0, 0, 0, 0], nrbin=30, nzbin=150))
    assert set(s.species) == {0, 1}


# ------------------------------------------------------------------ anisotropy from parameter files
@pytest.mark.parametrize("fmt,text", [
    ("toml", '[halo]\nm200 = 100.0\n[nbody]\nncomp = 1\nnpart_comp = [4000]\nfmass_comp = [1.0]\n'
             'form_comp = ["nfw"]\nrscale_comp = [1.0]\nrtrunc_comp = [1.5]\nistaper = true\nbeta_comp = [-0.25]\n'),
    ("nml", "&HALO\nm200=100.\nc200=10.\n/\n&NBODY\nncomp=1\nnpart_comp(1)=4000\nfmass_comp(1)=1.\n"
            "form_comp(1)='nfw'\nrscale_comp(1)=1.\nrtrunc_comp(1)=1.5\nraniso_comp(1)=2.0\nistaper=.true.\n/\n"),
])
def test_anisotropy_through_the_cli(tmp_path, fmt, text):
    path = tmp_path / ("p.toml" if fmt == "toml" else "p.txt")
    path.write_text(text)
    main(["nbody", str(path), "--outdir", str(tmp_path)])
    out = tmp_path / "nbody.gdt.hdf5"
    check_ic(out, {1})
    with h5py.File(out) as f:
        x, v = f["PartType1/Coordinates"][:].astype(float), f["PartType1/Velocities"][:].astype(float)
    r = np.linalg.norm(x, axis=1)
    vr = (v * x).sum(1) / r
    outer = r > np.median(r)
    beta = 1 - 0.5 * np.mean((v ** 2).sum(1)[outer] - vr[outer] ** 2) / np.mean(vr[outer] ** 2)
    assert (beta < -0.1) if fmt == "toml" else (beta > 0.2)       # tangential vs radially biased outskirts
