"""Stitching several systems into one initial-condition file (port + extension of merge_comps).

Each system is first rotated (``euler`` angles, degrees, z-x-z convention:
rotate by ``phi`` about z, then ``theta`` about x, then ``psi`` about z),
then translated (``position``) and boosted (``velocity``).  Systems may be
tidied beforehand with the old ``merge_comps`` options: ``taper_profile``,
``reduce_stellar_mass`` and ``excise_types``.  The stitched species are
concatenated type by type; per-particle masses are kept as a ``Masses`` block
whenever the masses of a type differ.  The AREPO background grid, if wanted,
is built *once* on the combined gas so the cells of the different systems
share one tree and one box.

``orbit_two_body`` supplies positions/velocities for a binary merger from an
initial separation, pericentre and eccentricity.
"""
from pathlib import Path

import numpy as np

from scipy.spatial import cKDTree

from . import backgroundgrid
from .constants import GAMMA, GGDT
from .io import write_hdf5
from .system import System


def rotation_matrix(euler_deg):
    phi, theta, psi = np.radians(euler_deg)

    def rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def rx(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    return rz(phi) @ rx(theta) @ rz(psi)


def place(system: System, position=(0, 0, 0), velocity=(0, 0, 0), euler=(0, 0, 0)) -> System:
    """Return a rotated, translated and boosted copy."""
    out = system.copy()
    R = rotation_matrix(euler)
    for s in out.species.values():
        s["pos"] = s["pos"].astype(np.float64) @ R.T + np.asarray(position, dtype=float)
        s["vel"] = s["vel"].astype(np.float64) @ R.T + np.asarray(velocity, dtype=float)
    out.meta.update(position=list(map(float, position)), velocity=list(map(float, velocity)),
                    euler=list(map(float, euler)))
    return out


# ---- the merge_comps tidying options -----------------------------------------
def taper_profile(system: System, mass_fraction, rng) -> System:
    """merge_comps -taper_profile: find the radius r_t enclosing ``mass_fraction`` of the
    system's mass, then keep each particle with probability exp(-(r-r_t)/r_t) (1 inside r_t).
    Particle masses are not rescaled, as in the Fortran."""
    out = system.copy()
    pos = np.concatenate([s["pos"] for s in out.species.values()]).astype(float)
    mass = np.concatenate([np.broadcast_to(s["mass"], (len(s["pos"]),)) for s in out.species.values()])
    r = np.linalg.norm(pos, axis=1)
    order = np.argsort(r)
    cum = np.cumsum(mass[order])
    rt = r[order][max(0, int(np.searchsorted(cum, mass_fraction * mass.sum())) - 1)]
    for t, s in out.species.items():
        rr = np.linalg.norm(s["pos"], axis=1)
        keep = rng.random(len(rr)) < np.exp(-(rr - rt) / rt)
        _subset(s, keep)
    out.meta["taper_radius"] = float(rt)
    return out


def reduce_stellar_mass(system: System, ptype, fraction) -> System:
    out = system.copy()
    if fraction == 0:
        out.species.pop(ptype, None)
    else:
        out.species[ptype]["mass"] = out.species[ptype]["mass"] * fraction
    return out


def excise_types(system: System, ptypes) -> System:
    out = system.copy()
    for t in ptypes:
        out.species.pop(t, None)
    return out


def _subset(s, keep):
    for k in ("pos", "vel", "ug"):
        if k in s:
            s[k] = s[k][keep]
    if np.ndim(s["mass"]):
        s["mass"] = s["mass"][keep]


# ---- orbit helper ---------------------------------------------------------------
def orbit_two_body(m1, m2, r_init, r_peri, ecc=1.0):
    """Initial conditions of a two-body orbit in the x-y plane, centre-of-mass frame.

    The pair starts incoming at separation ``r_init`` on a conic with
    pericentre ``r_peri`` and eccentricity ``ecc`` (1 = parabolic).
    Returns (pos1, vel1, pos2, vel2).
    """
    if r_init <= r_peri:
        raise ValueError("r_init must exceed the pericentre distance")
    mu = GGDT * (m1 + m2)
    p = r_peri * (1 + ecc)
    cos_t = (p / r_init - 1.0) / ecc
    if abs(cos_t) > 1:
        raise ValueError("orbit never reaches r_init (bound orbit with apocentre < r_init)")
    theta = -np.arccos(cos_t)                     # incoming branch
    vr = np.sqrt(mu / p) * ecc * np.sin(theta)
    vt = np.sqrt(mu / p) * (1 + ecc * np.cos(theta))
    c, s = np.cos(theta), np.sin(theta)
    rel_pos = r_init * np.array([c, s, 0.0])
    rel_vel = np.array([vr * c - vt * s, vr * s + vt * c, 0.0])
    M = m1 + m2
    return (-m2 / M * rel_pos, -m2 / M * rel_vel, m1 / M * rel_pos, m1 / M * rel_vel)


# ---- gas interfaces between systems ------------------------------------------------
class GasField:
    """Particle-based gas pressure estimate P = (gamma-1) rho u at arbitrary points.

    rho is the mass of the k nearest gas particles over the volume of the sphere
    reaching the k-th; u their mass-weighted mean.  The same estimator is used for
    a system's own particles and for queries from other systems, so the two are
    directly comparable (and it needs nothing but the particles, so it also works
    for systems loaded from files).
    """

    def __init__(self, pos, mass, u, k=32):
        self.k = int(min(k, len(pos)))
        self.tree = cKDTree(pos)
        self.mass = np.broadcast_to(np.asarray(mass, dtype=float), (len(pos),))
        self.u = np.asarray(u, dtype=float)

    def pressure(self, x):
        d, idx = self.tree.query(np.atleast_2d(x), k=self.k, workers=-1)
        d = d.reshape(len(np.atleast_2d(x)), -1)
        idx = idx.reshape(d.shape)
        r = np.maximum(d[:, -1], 1e-30)
        m = self.mass[idx].sum(1)
        rho = m / (4.0 / 3.0 * np.pi * r ** 3)
        u = (self.mass[idx] * self.u[idx]).sum(1) / m
        return (GAMMA - 1.0) * rho * u


def resolve_gas_overlap(placed, k=32, verbose=True):
    """Pressure-balance the interfaces where the gas of different systems overlaps.

    At each gas particle of system i, compare its own pressure P_i with the
    pressure P_j that every other system j would have at that point; the
    particle is removed if any P_j > P_i.  So in a region where two gas
    distributions overlap, the higher-pressure one is kept and the other is
    cut away, instead of the two being superposed (which would double the
    mass and put the interface out of pressure balance).  All pressures are
    taken from the *original* particle sets, so the rule is symmetric.
    Returns (new systems, stats).
    """
    log = print if verbose else (lambda *a, **k: None)
    gas = [s for s in placed if 0 in s.species and len(s.species[0]["pos"]) > 1]
    if len(gas) < 2:
        return placed, {}
    fields = {s.name: GasField(s.species[0]["pos"], s.species[0]["mass"], s.species[0]["ug"], k)
              for s in gas}
    out, stats = [], {}
    for s in placed:
        if s.name not in fields:
            out.append(s)
            continue
        g = s.species[0]
        p_own = fields[s.name].pressure(g["pos"])
        lose = np.zeros(len(p_own), dtype=bool)
        for other, f in fields.items():
            if other != s.name:
                lose |= f.pressure(g["pos"]) > p_own
        t = s.copy()
        m_before = float(np.sum(np.broadcast_to(g["mass"], (len(lose),))))
        _subset(t.species[0], ~lose)
        m_lost = m_before - float(np.sum(np.broadcast_to(t.species[0]["mass"], (len(t.species[0]["pos"]),))))
        stats[s.name] = dict(removed=int(lose.sum()), of=int(len(lose)), mass_removed=m_lost,
                             mass_fraction=m_lost / m_before)
        t.meta["gas_overlap"] = stats[s.name]
        out.append(t)
        log(f"[{s.name}] gas overlap: removed {int(lose.sum())} of {len(lose)} particles "
            f"({100 * m_lost / m_before:.2f}% of its gas mass)")
    return out, stats


# ---- stitching ------------------------------------------------------------------
def stitch(systems, type_maps=None, verbose=True):
    """Concatenate placed systems into one species dict.

    ``type_maps``: optional per-system dict {old_ptype: new_ptype} so e.g. two
    dark halos can be kept as separate GADGET types.
    """
    log = print if verbose else (lambda *a, **k: None)
    type_maps = type_maps or [None] * len(systems)
    parts = {}
    for sysm, tmap in zip(systems, type_maps):
        for t, s in sysm.species.items():
            nt = (tmap or {}).get(t, t)
            n = len(s["pos"])
            parts.setdefault(nt, []).append(dict(
                pos=s["pos"].astype(np.float64), vel=s["vel"].astype(np.float64),
                mass=np.broadcast_to(np.asarray(s["mass"], dtype=np.float64), (n,)),
                ug=s.get("ug")))
    merged = {}
    for t, lst in sorted(parts.items()):
        has_u = [p["ug"] is not None for p in lst]
        if any(has_u) and not all(has_u):
            raise ValueError(f"PartType{t}: only some contributing systems carry InternalEnergy")
        s = dict(pos=np.concatenate([p["pos"] for p in lst]),
                 vel=np.concatenate([p["vel"] for p in lst]),
                 mass=np.concatenate([p["mass"] for p in lst]))
        if lst[0]["ug"] is not None:
            s["ug"] = np.concatenate([p["ug"] for p in lst])
        if np.all(s["mass"] == s["mass"][0]):
            s["mass"] = float(s["mass"][0])
        merged[t] = s
        log(f"PartType{t}: {len(s['pos'])} particles from {len(lst)} system(s)")
    return merged


def finalise(species, background_grid=False, boxsize_factor=5.0, nngb=64, verbose=True, suggest_out=None):
    """Optionally convert the gas to an AREPO background-grid mesh.  Returns (species, boxsize).

    If ``suggest_out`` is a dict it receives the AREPO parameter suggestions (BoxSize, ReferenceGasPartMass, ...)."""
    if not background_grid:
        return species, 0.0
    if 0 not in species:
        raise ValueError("background grid requested but there is no gas (PartType0)")
    g = species[0]
    mass = np.broadcast_to(g["mass"], (len(g["pos"]),))
    bgo = backgroundgrid.add_background_grid(g["pos"], g["vel"], mass, g["ug"],
                                             boxsize_factor=boxsize_factor, nngb=nngb, verbose=verbose)
    out = {t: dict(s, pos=np.asarray(s["pos"], dtype=float) + bgo["offset"]) for t, s in species.items() if t != 0}
    out[0] = dict(pos=bgo["pos"], vel=bgo["vel"], mass=bgo["mass"], ug=bgo["u"])
    if verbose:
        print("AREPO suggestions:", bgo["suggest"])
    if suggest_out is not None:
        suggest_out.update(bgo["suggest"])
    return dict(sorted(out.items())), bgo["boxsize"]


def write_arepo_suggestions(ic_path, suggest):
    """Write the AREPO parameter suggestions next to the IC file as ``<ic_path>.arepo.json``."""
    import json
    if suggest:
        Path(str(ic_path) + ".arepo.json").write_text(json.dumps(suggest, indent=2))


def write_merged(path, species, boxsize=0.0, double=None, verbose=True):
    write_hdf5(path, species, boxsize=boxsize, double=bool(boxsize) if double is None else double,
               verbose=verbose)


def _extent(system, q=0.99):
    """Radius (about the system's own gas/particle centre) enclosing a fraction q of the gas
    (or of all particles if there is no gas)."""
    s = system.species.get(0) or next(iter(system.species.values()))
    c = np.asarray(system.meta.get("position", [0, 0, 0]), dtype=float)
    return float(np.quantile(np.linalg.norm(s["pos"] - c, axis=1), q))


def _warn_overlap(placed, verbose=True):
    """Systems are superposed, not re-equilibrated: flag pairs whose gas extents overlap."""
    if not verbose:
        return
    ext = [_extent(s) for s in placed]
    for i in range(len(placed)):
        for j in range(i + 1, len(placed)):
            sep = np.linalg.norm(np.asarray(placed[i].meta["position"]) - np.asarray(placed[j].meta["position"]))
            if sep < ext[i] + ext[j]:
                print(f"warning: '{placed[i].name}' and '{placed[j].name}' start {sep:.0f} kpc apart but "
                      f"their gas extends to {ext[i]:.0f} + {ext[j]:.0f} kpc; the overlapping gas is simply "
                      f"superposed (not in equilibrium). Increase r_init or lower rtrunc to avoid this.")


def merge_from_config(cfg_path, outdir=".", seed=None, verbose=True):
    """Build/load every system in a TOML file, place them and write one IC file.

    Example::

        output = "merger.hdf5"
        seed = 7
        [gas_overlap]                # optional; default mode = "pressure"
        mode = "pressure"            # or "none" (plain superposition, warns if extents overlap)
        k = 32                       # neighbours for the pressure estimate
        [background_grid]            # optional
        enabled = true
        boxsize_factor = 5
        nngb = 64
        [orbit]                      # optional: sets position/velocity of two systems
        systems = [0, 1]
        r_init = 300.0
        r_peri = 20.0
        ecc = 1.0
        [[system]]
        name = "gal1"
        params = "mw.txt"            # build from a parameter file ...
        euler = [0, 30, 0]
        taper_mass = 0.99            # optional merge_comps-style tidying
        [system.options]             # optional: kwargs of build_system
        halo_core = 0.1
        [system.overrides]           # optional: Params fields
        ndisc = 100000
        [[system]]
        name = "gal2"
        file = "other.hdf5"          # ... or load an existing file
        position = [300, 0, 0]
        velocity = [-50, 0, 0]
        type_map = { 1 = 2 }         # keep this system's halo as PartType2
    """
    import tomllib
    log = print if verbose else (lambda *a, **k: None)
    cfg_path = Path(cfg_path)
    with open(cfg_path, "rb") as fh:
        cfg = tomllib.load(fh)
    base = cfg_path.parent
    seed = cfg.get("seed", 1) if seed is None else seed
    rng = np.random.default_rng(seed)
    if not cfg.get("system"):
        raise ValueError("config has no [[system]] entries")

    from .system import build_system, read_system
    built = []
    for i, sc in enumerate(cfg["system"]):
        name = sc.get("name", f"system{i}")
        if "file" in sc:
            sysm = read_system(base / sc["file"], name)
        elif "params" in sc:
            sysm = build_system(base / sc["params"], name=name, seed=seed + 1000 * i + 1,
                                overrides=sc.get("overrides"), verbose=verbose, **sc.get("options", {}))
        else:
            raise ValueError(f"system {name!r} needs either 'params' or 'file'")
        if "taper_mass" in sc:
            sysm = taper_profile(sysm, sc["taper_mass"], rng)
        if "reduce_stellar_mass" in sc:
            r = sc["reduce_stellar_mass"]
            sysm = reduce_stellar_mass(sysm, int(r["type"]), float(r["fraction"]))
        if "excise_types" in sc:
            sysm = excise_types(sysm, [int(t) for t in sc["excise_types"]])
        log(f"[{name}] {sysm.counts()}, mass {sysm.mass():.5g}")
        built.append(sysm)

    pos = [np.array(sc.get("position", [0, 0, 0]), dtype=float) for sc in cfg["system"]]
    vel = [np.array(sc.get("velocity", [0, 0, 0]), dtype=float) for sc in cfg["system"]]
    if "orbit" in cfg:
        o = cfg["orbit"]
        i, j = o["systems"]
        m1 = o.get("m1", built[i].meta.get("m200", built[i].mass()))
        m2 = o.get("m2", built[j].meta.get("m200", built[j].mass()))
        pos[i], vel[i], pos[j], vel[j] = orbit_two_body(m1, m2, o["r_init"], o["r_peri"], o.get("ecc", 1.0))
        log(f"orbit: M1={m1:.4g} M2={m2:.4g}, r_init={o['r_init']}, r_peri={o['r_peri']}, e={o.get('ecc', 1.0)}")

    placed = [place(s, pos[k], vel[k], sc.get("euler", (0, 0, 0)))
              for k, (s, sc) in enumerate(zip(built, cfg["system"]))]
    ov = cfg.get("gas_overlap", {})
    mode = ov.get("mode", "pressure")
    if mode == "pressure":
        placed, _ = resolve_gas_overlap(placed, ov.get("k", 32), verbose)
    elif mode == "none":
        _warn_overlap(placed, verbose)
    else:
        raise ValueError(f"gas_overlap.mode must be 'pressure' or 'none', got {mode!r}")
    maps = [{int(a): int(b) for a, b in sc.get("type_map", {}).items()} or None for sc in cfg["system"]]
    species = stitch(placed, maps, verbose=verbose)
    bg = cfg.get("background_grid", {})
    suggest = {}
    species, boxsize = finalise(species, bg.get("enabled", False), bg.get("boxsize_factor", 5.0),
                                bg.get("nngb", 64), verbose=verbose, suggest_out=suggest)
    out = Path(outdir) / cfg.get("output", "merged.hdf5")
    out.parent.mkdir(parents=True, exist_ok=True)
    write_merged(out, species, boxsize, verbose=verbose)
    write_arepo_suggestions(out, suggest)
    return out
