#!/usr/bin/env python3
"""Set up a run directory for generating ICs (galaxy-components) and running Arepo on a SLURM cluster.

    python hpc/setup_run.py --galaxy hpc/cases/hot_halo/galaxy.toml --run hpc/cases/hot_halo/run.toml \\
        --cluster hpc/cluster.env --run-dir $SCRATCH/hot_halo [--scale 0.01] [--submit]

It writes into RUN_DIR:  galaxy.toml (the model), make_ics.sbatch, arepo.sbatch, submit.sh, run_settings.json,
and logs/.  The IC job builds the ICs with a background grid and then renders arepo_param.txt from the IC file's
header (BoxSize, gas particle mass), so the Arepo parameter file always matches the ICs it will read.  ``--submit``
runs submit.sh (the Arepo job depends on the IC job).  See hpc/README.md.

Only the standard library plus the package itself (for ``--scale``) and h5py (for the IC header) are needed.
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
UNIT_TIME_GYR = 3.085678e21 / 1e5 / 3.15576e16        # kpc / (km/s) in Gyr
CLUSTER_DEFAULTS = dict(
    ACCOUNT="", PARTITION="", MODULES="", ENV_SETUP="", LAUNCHER="srun", AREPO_BIN="",
    IC_CPUS="8", IC_MEM="64G", IC_TIME="04:00:00",
    AREPO_NODES="2", AREPO_NTASKS_PER_NODE="32", AREPO_MEM_PER_CPU="", AREPO_TIME="24:00:00",
)
RUN_DEFAULTS = dict(
    run=dict(name="run", snapshot="", t_end_gyr=5.0, n_snapshots=50, max_timestep_myr=5.0),
    ics=dict(background_grid=True, boxsize_factor=5.0, nngb=64, seed=19721),
    arepo=dict(softening_gas=0.1, softening_dm=0.5, softening_stars=0.1, param_template=""),
)


def read_env(path):
    """KEY=VALUE lines (shell style, optional quotes); a missing file gives the defaults."""
    out = dict(CLUSTER_DEFAULTS)
    if path:
        for line in Path(path).read_text().splitlines():
            m = re.match(r"\s*(?:export\s+)?([A-Z_][A-Z0-9_]*)=(.*)$", line)
            if m and not line.lstrip().startswith("#"):
                v = m.group(2).split(" #")[0].strip()
                out[m.group(1)] = v[1:-1] if len(v) > 1 and v[0] == v[-1] and v[0] in "\"'" else v
    return out


def read_run(path):
    out = {k: dict(v) for k, v in RUN_DEFAULTS.items()}
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    for table, vals in data.items():
        if table not in out:
            raise SystemExit(f"{path}: unknown table [{table}]; valid: {', '.join(out)}")
        for k, v in vals.items():
            if k not in out[table]:
                raise SystemExit(f"{path}: unknown field '{k}' in [{table}]; valid: {', '.join(out[table])}")
            out[table][k] = v
    return out


def render(text, values):
    """Substitute @KEY@; an #SBATCH line whose value is empty is dropped, any other unresolved key is an error."""
    lines = []
    for line in text.splitlines():
        keys = re.findall(r"@([A-Z0-9_]+)@", line)
        if line.startswith("#SBATCH") and any(values.get(k, None) in ("", None) for k in keys):
            continue
        for k in keys:
            if k not in values:
                raise SystemExit(f"template placeholder @{k}@ has no value")
            line = line.replace(f"@{k}@", str(values[k]))
        lines.append(line)
    return "\n".join(lines) + "\n"


def module_lines(cl):
    return "\n".join(f"module load {m}" for m in cl["MODULES"].split()) or "# (no modules requested)"


def ic_command(run, snapshot):
    ics = run["ics"]
    bg = f"--background-grid --boxsize-factor {ics['boxsize_factor']} --nngb {ics['nngb']}" if ics["background_grid"] else ""
    if not snapshot:
        return f"galaxy-components galaxy galaxy.toml --outdir ICs --name ic.hdf5 --seed {ics['seed']} {bg}".strip()
    return f"galaxy-components merge merge.toml --outdir ICs --seed {ics['seed']}"


def merge_toml(run, snapshot):
    ics = run["ics"]
    return f'''# gas (disc + hot halo, from galaxy.toml) added to an existing collisionless snapshot, e.g. from GalIC
output = "ic.hdf5"
seed = {ics['seed']}

[background_grid]
enabled = {str(bool(ics['background_grid'])).lower()}
boxsize_factor = {ics['boxsize_factor']}
nngb = {ics['nngb']}

[[system]]
name = "gas"
params = "galaxy.toml"
[system.options]
nbody_components = false      # the dark matter, bulge and stars come from the snapshot

[[system]]
name = "snapshot"
file = "{snapshot}"
'''


def scaled_galaxy(src, dst, scale):
    """Copy galaxy.toml, multiplying the particle numbers by ``scale`` (for quick test runs)."""
    if scale == 1.0:
        shutil.copy(src, dst)
        return
    try:
        from galaxy_components.params import read_params, to_toml
    except ImportError:                       # running from a checkout without installing the package
        sys.path.insert(0, str(HERE.parent))
        from galaxy_components.params import read_params, to_toml
    p = read_params(src)
    p.ndisc, p.nhalo = int(p.ndisc * scale), int(p.nhalo * scale)
    p.npart_comp = [int(n * scale) for n in p.npart_comp]
    Path(dst).write_text(to_toml(p, header=f"{Path(src).name} with particle numbers x {scale:g}"))


def setup(args):
    run_dir = Path(args.run_dir).resolve()
    cl, run = read_env(args.cluster), read_run(args.run)
    snapshot = str(Path(run["run"]["snapshot"]).expanduser().resolve()) if run["run"]["snapshot"] else ""
    if run["run"]["snapshot"] and not Path(snapshot).exists() and not args.allow_missing_snapshot:
        raise SystemExit(f"snapshot {snapshot} does not exist (set [run] snapshot in {args.run})")
    try:
        from galaxy_components.params import read_params
    except ImportError:
        sys.path.insert(0, str(HERE.parent))
        from galaxy_components.params import read_params
    p = read_params(args.galaxy)
    if not snapshot and not any(n > 0 for n in p.npart_comp[:p.ncomp]):
        raise SystemExit(f"{args.galaxy} has no N-body particles (npart_comp), so there would be no dark matter: "
                         f"set [run] snapshot in {args.run} to a collisionless snapshot, or give it N-body components")
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    scaled_galaxy(args.galaxy, run_dir / "galaxy.toml", args.scale)
    if snapshot:
        (run_dir / "merge.toml").write_text(merge_toml(run, snapshot))
    values = dict(cl, NAME=run["run"]["name"], RUN_DIR=run_dir, SETUP_PY=Path(__file__).resolve(),
                  MODULE_LINES=module_lines(cl), IC_COMMAND=ic_command(run, snapshot),
                  ENV_LINE=cl["ENV_SETUP"] or "# (no environment setup)")
    for tmpl, out in (("make_ics.sbatch.in", "make_ics.sbatch"), ("arepo.sbatch.in", "arepo.sbatch")):
        (run_dir / out).write_text(render((HERE / "templates" / tmpl).read_text(), values))
    (run_dir / "submit.sh").write_text(render((HERE / "templates" / "submit.sh.in").read_text(), values))
    (run_dir / "submit.sh").chmod(0o755)
    (run_dir / "run_settings.json").write_text(json.dumps(dict(cluster=cl, run=run, snapshot=snapshot), indent=2))
    print(f"set up {run_dir}")
    if not cl["AREPO_BIN"]:
        print("note: AREPO_BIN is not set in the cluster file; arepo.sbatch will refuse to run until it is")
    if args.submit:
        subprocess.run(["bash", str(run_dir / "submit.sh")], check=True)
    else:
        print(f"next: cd {run_dir} && ./submit.sh      (or sbatch make_ics.sbatch, then sbatch arepo.sbatch)")


# ------------------------------------------------------------------ Arepo parameter file from the ICs
def ic_header(path):
    """BoxSize, particle counts and the reference gas particle mass of an IC file.

    With a background grid the gas is cells of unequal mass, so the mean in the file is not the particle mass the
    cells should be refined towards; galaxy-components writes that mass (and the box) to ``<ic>.arepo.json``."""
    import h5py
    import numpy as np
    with h5py.File(path, "r") as f:
        h = f["Header"].attrs
        n = np.asarray(h["NumPart_ThisFile"], dtype=int)
        mt = np.asarray(h["MassTable"], dtype=float)
        gas_mass = float(mt[0]) if mt[0] > 0 else (float(f["PartType0"]["Masses"][:].mean()) if n[0] else 0.0)
        out = dict(boxsize=float(h["BoxSize"]), npart=n.tolist(), gas_mass=gas_mass)
    side = Path(str(path) + ".arepo.json")
    if side.exists():
        sug = json.loads(side.read_text())
        out["gas_mass"] = float(sug.get("ReferenceGasPartMass", gas_mass))
    return out


def render_arepo(run_dir):
    run_dir = Path(run_dir).resolve()
    st = json.loads((run_dir / "run_settings.json").read_text())
    cl, run = st["cluster"], st["run"]
    hdr = ic_header(run_dir / "ICs" / "ic.hdf5")
    if run["ics"]["background_grid"] and hdr["boxsize"] <= 0:
        raise SystemExit("the ICs have BoxSize 0 although a background grid was requested")
    t_end = run["run"]["t_end_gyr"] / UNIT_TIME_GYR
    walltime = sum(int(x) * m for x, m in zip(cl["AREPO_TIME"].replace("-", ":").split(":")[::-1], (1, 60, 3600, 86400)))
    soft = run["arepo"]
    tmpl = Path(soft["param_template"]) if soft["param_template"] else HERE / "templates" / "arepo_param.txt.in"
    values = dict(
        IC_FILE=str(run_dir / "ICs" / "ic"), OUTPUT_DIR=str(run_dir / "output"), BOXSIZE=f"{hdr['boxsize']:.8g}",
        TIME_MAX=f"{t_end:.6g}", TIME_BET_SNAP=f"{t_end / run['run']['n_snapshots']:.6g}",
        MAX_DT=f"{run['run']['max_timestep_myr'] * 1e-3 / UNIT_TIME_GYR:.6g}",
        TIMELIMIT_CPU=max(walltime - 900, 600), REF_GAS_MASS=f"{hdr['gas_mass']:.8g}",
        SOFT_GAS=soft["softening_gas"], SOFT_DM=soft["softening_dm"], SOFT_STARS=soft["softening_stars"],
        NWRITE=min(8, int(cl["AREPO_NODES"]) * int(cl["AREPO_NTASKS_PER_NODE"])))
    (run_dir / "arepo_param.txt").write_text(render(tmpl.read_text(), values))
    print(f"wrote {run_dir / 'arepo_param.txt'}: BoxSize {hdr['boxsize']:.6g} kpc, particles {hdr['npart']}, "
          f"gas particle mass {hdr['gas_mass']:.4g}, TimeMax {t_end:.4g} ({run['run']['t_end_gyr']:g} Gyr)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxy", help="galaxy-components parameter file (TOML)")
    ap.add_argument("--run", help="run settings (TOML): snapshot, duration, softenings")
    ap.add_argument("--cluster", help="cluster settings (KEY=VALUE): account, modules, Arepo binary, resources")
    ap.add_argument("--run-dir")
    ap.add_argument("--scale", type=float, default=1.0, help="multiply all particle numbers (quick test runs)")
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--allow-missing-snapshot", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--render-arepo", metavar="RUN_DIR", help="(used by the IC job) write arepo_param.txt from the ICs")
    args = ap.parse_args(argv)
    if args.render_arepo:
        return render_arepo(args.render_arepo)
    if not (args.galaxy and args.run and args.run_dir):
        ap.error("--galaxy, --run and --run-dir are required")
    setup(args)


if __name__ == "__main__":
    main()
