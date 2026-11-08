"""Command-line entry points:  python -m galaxy_components {disc,nbody} parameters.txt"""
import argparse
from pathlib import Path

import numpy as np

from . import backgroundgrid, compound, gas_disc, gas_halo, nbody
from .constants import PI, GGDT
from .io import write_hdf5
from .model import build_model
from .params import output_name, read_params, to_toml


def run_disc(args):
    p = read_params(args.params)
    model = build_model(p)
    disc = gas_disc.build_gas_disc(p, model)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    _write_disc_tables(out, p, model, disc)
    rng = np.random.default_rng(args.seed)
    pos, vel, ug, mpart = gas_disc.sample_particles(p, disc, p.ndisc, rng)
    species = {0: dict(pos=pos, vel=vel, ug=ug, mass=mpart)}
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)
    write_hdf5(out / output_name(p.disc_file, p.snap_format), species)


def _write_disc_tables(out, p, model, disc):
    r = disc.rtab
    hdr = (f"rdisc [kpc]: {model.rdisc:.8f}\nmdisc [1e10 Msol]: {model.mdisc:.8f}\n"
           f"mstar [1e10 Msol]: {model.mstar:.8f}\nsigma0 [1e10 Msol/kpc^2]: {disc.sigma0:.8f}\n"
           "Radius[kpc] Vc_bg[km/s] Vdisc[km/s] Vp[km/s] Vrot[km/s]")
    v = disc.vcomp
    np.savetxt(out / "vrot_components.txt",
               np.column_stack((r, np.sqrt(v["background"]), np.sqrt(v["gas"]),
                                np.sqrt(np.abs(v["pressure"])), disc.vrot)), header=hdr)
    dv = np.gradient(disc.vrot, r)
    kappa = np.sqrt(np.maximum(2 * (disc.vrot / r) * (dv + disc.vrot / r), 0.0))
    sig = disc.sigma0 * np.exp(-r / model.rdisc)
    q = disc.cs * kappa / PI / GGDT / sig
    np.savetxt(out / "disc_parameters.txt",
               np.column_stack((r, sig, disc.vrot, q, sig / 2 / disc.rho0)),
               header=hdr.split("\n")[0] + "\nRadius Sigma(R) Vrot Q Sigma/2rho0")


def run_nbody(args):
    p = read_params(args.params)
    nb = nbody.build_nbody(p, legacy=not args.fix_legacy, egrid="log" if args.log_egrid else None)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rscale0 = nb.rscale_comp[0]
    np.savetxt(out / "distribution_function.txt",
               np.column_stack((nb.etab, rscale0 ** 3 * nb.ftab.T)),
               header="Binding Energy E [g*mtot/rscale]  F(E)")
    rng = np.random.default_rng(args.seed)
    comps = nbody.sample_nbody(p, nb, rng)
    s = nbody.summarise(nb, comps)
    print(f"Generated halo with {sum(len(c['pos']) for c in comps)} particles; {s['nvir']} inside r200")
    print("Centre of mass [kpc]:", s["com"], "\nCentre of mass velocity [km/s]:", s["vcom"])
    print("Velocity dispersion [km/s]:", s["sigma"], s["sigma_tot"])
    print("Radial velocity (dispersion) [km/s]:", s["vr"], s["sigma_r"])
    species = {i + 1: dict(pos=c["pos"], vel=c["vel"], mass=c["mpart"]) for i, c in enumerate(comps)}
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)
    write_hdf5(out / output_name(p.nbody_file, p.snap_format), species)


def _core(p, args):
    return args.halo_core if args.halo_core is not None else p.halo_core


def run_halo(args):
    p = read_params(args.params)
    from .system import build_halo
    halo, pg = build_halo(p, legacy_halo=args.legacy, halo_core=_core(p, args),
                          zero_edge=args.zero_edge or p.halo_zero_edge)
    rng = np.random.default_rng(args.seed)
    pos, vel, ug, mpart = gas_halo.sample_gas_halo(pg, halo, rng, velocity_noise=not args.cold)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    species = {0: dict(pos=pos, vel=vel, ug=ug, mass=mpart)}
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)
    write_hdf5(out / output_name(p.halo_file, p.snap_format), species)


def run_galaxy(args):
    """One system (disc + halo gas + N-body) in one file, optionally with an AREPO background grid."""
    from .merge import finalise, write_merged
    from .system import build_system
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    p = read_params(args.params)
    sysm = build_system(p, name="galaxy", seed=args.seed, naive_merge=args.naive_merge,
                        halo_core=_core(p, args) if args.halo_core is not None else None,
                        zero_edge=True if args.zero_edge else None,
                        legacy_halo=args.legacy_halo, cold_halo=args.cold_halo,
                        nbody_components=not args.no_nbody, verbose=True)
    species, boxsize = finalise(sysm.species, args.background_grid, args.boxsize_factor, args.nngb)
    write_merged(out / (args.name or "galaxy_ic.hdf5"), species, boxsize)


def run_system(args):
    """Alias of ``galaxy``: build a single system file that can later be stitched."""
    run_galaxy(args)


def run_convert(args):
    """Convert a parameter file (e.g. a legacy Fortran namelist) to TOML."""
    p = read_params(args.params)
    text = to_toml(p, header=f"Converted from {Path(args.params).name} by galaxy-components convert.\n"
                             "Only values that differ from the defaults are written.", all_fields=args.all)
    if args.output:
        Path(args.output).write_text(text)
        print(f"wrote {args.output}")
    else:
        print(text, end="")


def run_check(args):
    """Background model, gas-disc rotation curve and vertical resolution, without sampling particles."""
    from . import workflow
    p = read_params(args.params)
    if args.ndisc:
        p.ndisc = args.ndisc
    print(workflow.model_report(p).format())


def run_add_gas(args):
    """Sample the gas disc and add it to an existing (e.g. GalIC) snapshot."""
    from . import workflow
    p = read_params(args.params)
    if args.ndisc:
        p.ndisc = args.ndisc
    if args.sampling:
        p.sampling = args.sampling
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    workflow.add_gas_to_snapshot(p, args.snapshot, out, seed=args.seed, background_grid=args.background_grid,
                                 boxsize_factor=args.boxsize_factor, nngb=args.nngb, strict=args.strict)
    print(f"wrote {out}")


def run_merge(args):
    from .merge import merge_from_config
    merge_from_config(args.config, args.outdir, seed=args.seed)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="galaxy_components")
    sub = ap.add_subparsers(dest="cmd", required=True)
    cp = sub.add_parser("convert", help="convert a parameter file (e.g. a Fortran namelist) to TOML")
    cp.add_argument("params")
    cp.add_argument("-o", "--output", help="output .toml file (default: print to stdout)")
    cp.add_argument("--all", action="store_true", help="write every field, not only non-defaults")
    cp.set_defaults(fn=run_convert)
    kp = sub.add_parser("check", help="gas-disc model: rotation curve, scale heights and vertical resolution")
    kp.add_argument("params")
    kp.add_argument("--ndisc", type=int, help="particle number to assess (default: the parameter file's)")
    kp.set_defaults(fn=run_check)
    ap_ = sub.add_parser("add-gas", help="add the gas disc to an existing snapshot (e.g. from GalIC)")
    ap_.add_argument("params")
    ap_.add_argument("--snapshot", required=True, help="GADGET HDF5 file without gas")
    ap_.add_argument("--out", required=True, help="output IC file")
    ap_.add_argument("--ndisc", type=int)
    ap_.add_argument("--sampling", choices=("poisson", "sobol"))
    ap_.add_argument("--seed", type=int, default=19721)
    ap_.add_argument("--strict", action="store_true", help="fail if the snapshot disagrees with the analytic model")
    ap_.add_argument("--background-grid", action="store_true", help="AREPO background grid")
    ap_.add_argument("--boxsize-factor", type=float, default=5.0)
    ap_.add_argument("--nngb", type=int, default=64)
    ap_.set_defaults(fn=run_add_gas)
    mp = sub.add_parser("merge", help="stitch several systems (TOML config) into one IC file")
    mp.add_argument("config")
    mp.add_argument("--outdir", default=".")
    mp.add_argument("--seed", type=int, default=None)
    mp.set_defaults(fn=run_merge)
    for name, fn in (("disc", run_disc), ("nbody", run_nbody), ("halo", run_halo), ("galaxy", run_galaxy),
                     ("system", run_system)):
        sp = sub.add_parser(name)
        sp.add_argument("params")
        sp.add_argument("--outdir", default=".")
        sp.add_argument("--seed", type=int, default=19721)
        if name == "nbody":
            sp.add_argument("--fix-legacy", action="store_true",
                            help="correct Fortran edge-extrapolation quirk instead of reproducing it")
        if name in ("halo", "galaxy", "system"):
            sp.add_argument("--halo-core", type=float, default=None,
                            help="gas entropy floor / K&S entropy at r200 (overrides &MODEL halo_core)")
            sp.add_argument("--log-egrid", action="store_true",
                            help="Fortran energy grid (logarithmic in E) instead of one node per radial node")
            sp.add_argument("--zero-edge", action="store_true",
                            help="hydrostatic halo to rtrunc*r200 with zero pressure at the edge (Donnert 2013)")
        if name == "halo":
            sp.add_argument("--legacy", action="store_true",
                            help="Fortran halo potential (a=r200/c200, M200) instead of the disc's dark halo")
            sp.add_argument("--cold", action="store_true",
                            help="no Gaussian velocity noise (Fortran adds sigma=sqrt(2u))")
        if name in ("galaxy", "system"):
            sp.add_argument("--legacy-halo", action="store_true")
            sp.add_argument("--cold-halo", action="store_true")
            sp.add_argument("--naive-merge", action="store_true",
                            help="plain concatenation (merge_comps behaviour): no pressure-balance")
            sp.add_argument("--no-nbody", action="store_true")
            sp.add_argument("--background-grid", action="store_true", help="AREPO background grid")
            sp.add_argument("--boxsize-factor", type=float, default=5.0)
            sp.add_argument("--nngb", type=int, default=64)
            sp.add_argument("--name")
        sp.set_defaults(fn=fn)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
