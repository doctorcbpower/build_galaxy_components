"""The GalIC + gas-disc workflow: check the analytic model, build the gas disc, add it to a snapshot.

The usual use of this package: a collisionless galaxy (dark halo, bulge, stellar disc) already exists
-- typically from GalIC -- and the gas disc is added in the analytic potential of those components.
Three steps, each available from Python (these functions, used by the notebook) and from the command
line (``galaxy-components check`` / ``disc`` / ``add-gas``):

1. ``model_report``      the background model, gas-disc rotation curve and vertical resolution
2. ``compare_snapshot``  does the snapshot agree with the analytic model the gas disc assumes?
3. ``add_gas_to_snapshot``  sample the gas disc and write snapshot + gas as one IC file
"""
from dataclasses import dataclass, field

import numpy as np

from . import gas_disc
from .io import write_hdf5
from .merge import finalise, stitch, write_arepo_suggestions, write_merged
from .model import build_model
from .system import System, read_system

RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0)       # in units of the disc scale length
HALF_MASS_SECH2 = np.arctanh(0.5)             # half-mass |z| / z0 of a sech^2 sheet (= 0.549)


@dataclass
class Report:
    p: object
    model: object
    disc: object
    mpart: float
    rows: list
    summary: dict
    notes: list = field(default_factory=list)

    def format(self) -> str:
        m, d, p = self.model, self.disc, self.p
        s = self.summary
        gas = d.gas_mass
        lines = [
            f"M200 = {m.m200:.4g} (1e10 Msol), R200 = {m.r200:.1f} kpc, halo scale a = {m.halo_a:.2f} kpc",
            f"stars: disc {m.mstar:.4g} (R_d = {m.rdisc:.2f} kpc, z0/R_d = {p.zdisc:g}), bulge {m.mbulge:.4g}; "
            f"gas disc {gas:.4g} ({gas / max(m.mstar, 1e-30):.2f} of the stellar disc), cs = {d.cs:.2f} km/s",
            f"ndisc = {p.ndisc:,} ({p.sampling}), particle mass {self.mpart:.3e} (1e10 Msol)",
            "",
            "   R/R_d   R [kpc]   vc_bg   vc_gas   vrot   [km/s]   gas half-mass height [pc]   h/H",
        ]
        for r in self.rows:
            flag = "  <-- unresolved" if r["h_over_H"] > gas_disc.HRES_MAX else ""
            lines.append(f"  {r['R_over_Rd']:6.1f}  {r['R']:8.2f}  {r['vc_bg']:6.1f}  {r['vc_gas']:6.1f}  {r['vrot']:6.1f}"
                         f"      {1e3 * r['H']:14.1f}               {r['h_over_H']:5.1f}{flag}")
        lines.append("")
        if s["frac"] >= 0.05:
            lines.append(f"Vertical resolution: h/H > {gas_disc.HRES_MAX:g} for {100 * s['frac']:.0f}% of the gas mass; "
                         f"ndisc ~ {s['n90']:,} resolves 90% of it (h/H = {s['q90']:.1f} at the 90th mass percentile now).")
        else:
            lines.append(f"Vertical resolution: h/H <= {gas_disc.HRES_MAX:g} for {100 * (1 - s['frac']):.0f}% of the gas mass.")
        lines += self.notes
        return "\n".join(lines)


def model_report(p, ndisc=None) -> Report:
    """Build the background model and the gas disc for ``p`` (no particles) and summarise them."""
    model = build_model(p, verbose=False)
    if model.mstar <= 0 and model.mdisc <= 0:
        raise ValueError("no disc: set mdisc > 0")
    disc = gas_disc.build_gas_disc(p, model, verbose=False)
    n = ndisc or p.ndisc or 100000
    r0, rcut, m0, mkeep = gas_disc.disc_window(p, disc)
    mpart = (mkeep - m0) / n
    summary = gas_disc.resolution_summary(disc, mpart, n, r0, rcut)
    ratio = gas_disc.resolution(disc, mpart)["ratio"]
    Hz = np.array([np.interp(0.5, mz, disc.ztab) for mz in disc.mztab])
    rows = []
    for f in RADII:
        i = int(np.argmin(abs(disc.rtab - f * model.rdisc)))
        rows.append(dict(R_over_Rd=f, R=disc.rtab[i], vc_bg=float(np.sqrt(disc.vcomp["background"][i])),
                         vc_gas=float(np.sqrt(disc.vcomp["gas"][i])), vrot=float(disc.vrot[i]),
                         H=float(Hz[i]), h_over_H=float(ratio[i])))
    notes = []
    gas_mass = model.mdisc - model.mstar
    if gas_mass > 0 and p.halo_type == "hernquist":
        mh_galic = model.m200 - model.mbulge - model.mstar - model.mbh
        notes.append(f"GalIC convention: its halo has M200 - M_bulge - M_disc = {mh_galic:.4g}; here the gas disc is also taken "
                     f"out of M200, so the halo is {mh_galic - model.mhalo:.3g} ({100 * (mh_galic - model.mhalo) / mh_galic:.2f}%) lighter.")
    return Report(p, model, disc, mpart, rows, summary, notes)


# ------------------------------------------------------------------ snapshot checks
@dataclass
class Finding:
    name: str
    snapshot: float
    model: float
    ok: bool
    note: str = ""

    def format(self):
        return f"{'ok  ' if self.ok else 'WARN'}  {self.name:<28s} snapshot {self.snapshot:10.4g}   model {self.model:10.4g}   {self.note}"


def _mass(s):
    return float(np.sum(s["mass"])) if np.ndim(s["mass"]) else float(s["mass"]) * len(s["pos"])


def compare_snapshot(p, model, species) -> list:
    """Check a snapshot (GADGET types 1 halo, 2 stellar disc, 3 bulge) against the analytic model.

    The gas disc is built in the model's potential, so a mismatch in mass, scale length, thickness or
    spin direction means the gas will not be in equilibrium in the snapshot.
    """
    out = []

    def ratio(name, snap, mod, tol, note=""):
        out.append(Finding(name, snap, mod, abs(snap / mod - 1) <= tol if mod > 0 else snap == 0, note))

    if 1 in species:
        ratio("halo mass [1e10 Msol]", _mass(species[1]), model.m200 - model.mbulge - model.mstar - model.mbh, 0.05,
              "(GalIC: M200 - bulge - stellar disc)")
    if 3 in species:
        ratio("bulge mass [1e10 Msol]", _mass(species[3]), model.mbulge, 0.05)
    if 2 in species:
        d = species[2]
        ratio("stellar disc mass [1e10 Msol]", _mass(d), model.mstar, 0.05)
        R = np.hypot(d["pos"][:, 0], d["pos"][:, 1])
        w = np.broadcast_to(d["mass"], (len(R),))
        ratio("disc scale length [kpc]", float(np.average(R, weights=w) / 2), model.rdisc, 0.10, "(mean R / 2)")
        z0 = float(np.median(np.abs(d["pos"][:, 2])) / HALF_MASS_SECH2)
        want = p.zdisc * model.rdisc
        if want > 0:
            ratio("disc thickness z0 [kpc]", z0, want, 0.25, "(median |z| / 0.549; sech^2)")
        else:
            out.append(Finding("disc thickness z0 [kpc]", z0, 0.0, z0 < 0.02 * model.rdisc,
                               "zdisc = 0 in the model: razor-thin" if z0 >= 0.02 * model.rdisc else ""))
        lz = float(np.sum(w * (d["pos"][:, 0] * d["vel"][:, 1] - d["pos"][:, 1] * d["vel"][:, 0])))
        out.append(Finding("disc spin direction (Lz > 0)", lz, 1.0, lz > 0, "gas rotates counter-clockwise about +z"))
    ref = species[2] if 2 in species else None          # the gas disc sits at the origin: compare the stellar disc
    if ref is not None:
        m_ = np.broadcast_to(ref["mass"], (len(ref["pos"]),))
        com = (ref["pos"].astype(float) * m_[:, None]).sum(0) / m_.sum()
        out.append(Finding("stellar-disc centre offset [kpc]", float(np.linalg.norm(com)), 0.0,
                           np.linalg.norm(com) < 0.1 * model.rdisc, "the gas disc is centred on the origin"))
    return out


# ------------------------------------------------------------------ build and merge
def add_gas_to_snapshot(p, snapshot, out, seed=19721, background_grid=False, boxsize_factor=5.0, nngb=64,
                        strict=False, verbose=True) -> dict:
    """Sample the gas disc of ``p`` and write ``snapshot`` plus the gas as one IC file ``out``.

    ``snapshot`` is a GADGET HDF5 file (e.g. from GalIC) without gas.  With ``strict`` any failed
    ``compare_snapshot`` check is an error; otherwise they are returned (and printed if ``verbose``).
    """
    log = print if verbose else (lambda *a, **k: None)
    if p.ndisc <= 0:
        raise ValueError("ndisc must be > 0 to build a gas disc")
    snap = read_system(snapshot, "snapshot")
    if 0 in snap.species:
        raise ValueError(f"{snapshot} already contains gas (PartType0)")
    model = build_model(p, verbose=False)
    findings = compare_snapshot(p, model, snap.species)
    for f in findings:
        log(f.format())
    if strict and not all(f.ok for f in findings):
        raise ValueError("snapshot does not match the analytic model:\n" + "\n".join(f.format() for f in findings if not f.ok))
    disc = gas_disc.build_gas_disc(p, model, verbose=False)
    pos, vel, ug, mpart = gas_disc.sample_particles(p, disc, p.ndisc, np.random.default_rng(seed))
    gas = System("gas_disc", {0: dict(pos=pos, vel=vel, ug=ug, mass=mpart)})
    species = stitch([snap, gas], verbose=verbose)
    suggest = {}
    species, boxsize = finalise(species, background_grid, boxsize_factor, nngb, verbose=verbose, suggest_out=suggest)
    write_merged(out, species, boxsize, verbose=verbose)
    write_arepo_suggestions(out, suggest)
    return dict(out=out, findings=findings, counts={t: len(s["pos"]) for t, s in species.items()}, boxsize=boxsize)


def standin_snapshot(p, path, n_halo=40000, n_disc=20000, n_bulge=8000, seed=3):
    """Write a *stand-in* for a GalIC snapshot so the workflow can be run without one.

    A spherical Hernquist halo and bulge plus an exponential/sech^2 stellar disc with the model's
    own scale lengths and masses.  The velocities are crude (circular speed for the disc, half of it
    as the dispersion of the spheres), so this is **not** an equilibrium model: use real GalIC output.
    """
    rng = np.random.default_rng(seed)
    model = build_model(p, verbose=False)
    species = {}
    rg = np.logspace(-3, np.log10(2 * model.r200), 400)
    vgrid = np.sqrt(model.vc2_background(rg))
    vcirc = lambda r: np.interp(r, rg, vgrid)

    def hernquist(n, mass, a):
        rmax = model.r200
        qmax = (rmax / (rmax + a)) ** 2
        q = np.sqrt(rng.random(n) * qmax)
        r = a * q / (1 - q)
        u = rng.normal(size=(n, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        vc = vcirc(r)
        return dict(pos=u * r[:, None], vel=rng.normal(size=(n, 3)) * 0.5 * vc[:, None] / np.sqrt(3), mass=mass / n)

    mh = model.m200 - model.mbulge - model.mstar - model.mbh
    species[1] = hernquist(n_halo, mh, model.halo_a)
    if model.mbulge > 0:
        species[3] = hernquist(n_bulge, model.mbulge, p.fbulge * model.halo_a)
    R = rng.gamma(2.0, model.rdisc, n_disc)
    z0 = p.zdisc * model.rdisc
    z = z0 * np.arctanh(2 * rng.random(n_disc) - 1) if z0 > 0 else np.zeros(n_disc)
    ph = 2 * np.pi * rng.random(n_disc)
    vc = vcirc(R)
    species[2] = dict(pos=np.column_stack((R * np.cos(ph), R * np.sin(ph), z)),
                      vel=np.column_stack((-vc * np.sin(ph), vc * np.cos(ph), np.zeros(n_disc))), mass=model.mstar / n_disc)
    write_hdf5(path, species, verbose=False)
    return path
