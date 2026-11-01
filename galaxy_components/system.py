"""A self-contained galaxy ("system") built in memory, and a reader for existing files.

A System is a dict of particle species in GADGET type numbering,

    species[ptype] = dict(pos (N,3), vel (N,3), mass (scalar or (N,)), ug (N,) [gas only])

with positions/velocities in the system's own frame (centre at the origin, at
rest).  ``build_system`` assembles everything one parameter file describes --
gas disc, gas halo (compounded by pressure balance, see ``compound``) and the
spherical N-body components -- and ``read_system`` loads a system someone else
made (including files from the Fortran codes) so it can be stitched too.
"""
from dataclasses import dataclass, field, replace
from pathlib import Path

import h5py
import numpy as np

from . import compound, gas_disc, gas_halo, nbody
from .model import build_model
from .params import Params, read_params
from .potentials import ExponentialDisc


@dataclass
class System:
    name: str
    species: dict
    meta: dict = field(default_factory=dict)

    def copy(self):
        sp = {t: {k: np.array(v, copy=True) if isinstance(v, np.ndarray) else v for k, v in s.items()}
              for t, s in self.species.items()}
        return System(self.name, sp, dict(self.meta))

    def counts(self):
        return {t: len(s["pos"]) for t, s in self.species.items()}

    def mass(self, ptype=None):
        types = [ptype] if ptype is not None else list(self.species)
        return float(sum(np.sum(self.species[t]["mass"]) if np.ndim(self.species[t]["mass"])
                         else self.species[t]["mass"] * len(self.species[t]["pos"]) for t in types))

    def centre_of_mass(self):
        m = np.concatenate([np.broadcast_to(s["mass"], (len(s["pos"]),)) for s in self.species.values()])
        x = np.concatenate([s["pos"] for s in self.species.values()]).astype(float)
        return (m[:, None] * x).sum(0) / m.sum()


def build_halo(p, nb=None, legacy_halo=False, halo_core=None, zero_edge=None):
    """The gas halo for ``p`` and the Params it was built with.

    With ``istaper`` the N-body halo is truncated at ``rtrunc_comp * r200`` (the largest over the
    components), and the gas is truncated at the same radius.  The K&S temperature profile is the
    equilibrium of the *untruncated* potential; in the shallower tapered potential it would be too
    hot and the outskirts would expand.  So the gas is re-solved in the actual enclosed mass of the
    N-body model: K&S entropy shape (plus ``halo_core``), mass fbaryon*M200 inside r200, and P = 0
    at the truncation radius.  The pressure and temperature therefore fall to zero at the edge and
    the gas is in hydrostatic equilibrium in the tapered halo.  ``legacy_halo`` skips this, as the
    Fortran did.
    """
    halo, pg = _spherical_halo(p, nb, legacy_halo, halo_core, zero_edge)
    if p.halo_feels_disc and not legacy_halo:
        model = build_model(p, verbose=False)
        if model.mdisc > 0:
            halo = gas_halo.adapt_to_disc(halo, ExponentialDisc(model.mdisc, model.rdisc), verbose=False)
    return halo, pg


def _spherical_halo(p, nb, legacy_halo, halo_core, zero_edge):
    core = p.halo_core if halo_core is None else halo_core
    edge0 = p.halo_zero_edge if zero_edge is None else zero_edge
    if not (p.istaper and p.ncomp > 0) or legacy_halo:
        ks = gas_halo.build_gas_halo(p, consistent=not legacy_halo, verbose=False)
        return gas_halo.add_entropy_core(p, ks, core, verbose=False, zero_edge=edge0), p
    redge = max(p.rtrunc_comp[:p.ncomp])
    if redge <= 1.0:
        raise ValueError(f"the gas halo is truncated at the N-body taper radius, which must lie beyond "
                         f"r200 (max rtrunc_comp = {redge:g}); raise rtrunc_comp or set istaper = false")
    if p.rtrunc != Params().rtrunc and abs(p.rtrunc - redge) > 1e-9:
        raise ValueError(f"rtrunc = {p.rtrunc:g} conflicts with the N-body taper radius {redge:g} r200: "
                         "with istaper the gas halo ends where the N-body halo is truncated, so leave "
                         "rtrunc unset")
    pg = replace(p, rtrunc=redge)
    nb = nbody.build_nbody(p, verbose=False) if nb is None else nb
    lr, lm = np.log(nb.rtab), np.log(nb.mtab)
    ks = gas_halo.build_gas_halo(pg, consistent=True, verbose=False)
    halo = gas_halo.add_entropy_core(pg, ks, core, verbose=False, zero_edge=True,
                                     mdm=lambda r: np.exp(np.interp(np.log(r), lr, lm)))
    return halo, pg


def build_system(params, name="system", seed=1, overrides=None, naive_merge=False,
                 halo_core=None, zero_edge=None, legacy_halo=False, cold_halo=False, nbody_components=True,
                 verbose=False) -> System:
    """Build one galaxy.  ``params`` is a path or a Params; ``overrides`` a dict of Params fields."""
    p = read_params(params) if not isinstance(params, Params) else params
    for k, v in (overrides or {}).items():
        if not hasattr(p, k):
            raise ValueError(f"unknown parameter override {k!r}")
        setattr(p, k, v)
    rng = np.random.default_rng(seed)
    log = print if verbose else (lambda *a, **k: None)
    meta = {}
    species = {}

    have_disc = p.ndisc > 0 and p.mdisc != 0
    have_halo = p.nhalo > 0 and p.fbaryon > 0
    model = build_model(p, verbose=False)
    meta.update(m200=model.m200, r200=model.r200)

    gas_parts = []
    halo = disc = None
    if have_halo:
        halo, pg = build_halo(p, legacy_halo=legacy_halo, halo_core=halo_core, zero_edge=zero_edge)
        hpos, hvel, hug, hmp = gas_halo.sample_gas_halo(pg, halo, rng, velocity_noise=not cold_halo,
                                                        verbose=False)
    if have_disc:
        disc = gas_disc.build_gas_disc(p, model, verbose=False)
        sdisc = disc
        if have_halo and not naive_merge:
            sdisc, mfrac, _ = compound.truncate_disc(disc, halo)
            keep = compound.excise_halo(hpos, disc, halo)
            meta["compound"] = dict(compound.mass_summary(disc, sdisc, mfrac),
                                    halo_excised=int((~keep).sum()))
            log(f"[{name}] {meta['compound']}")
            hpos, hvel, hug = hpos[keep], hvel[keep], hug[keep]
        dpos, dvel, dug, dmp = gas_disc.sample_particles(p, sdisc, p.ndisc, rng)
        gas_parts.append((dpos, dvel, dug, np.full(len(dpos), dmp)))
    if have_halo:
        gas_parts.append((hpos, hvel, hug, np.full(len(hpos), hmp)))
    if gas_parts:
        species[0] = dict(pos=np.concatenate([g[0] for g in gas_parts]),
                          vel=np.concatenate([g[1] for g in gas_parts]),
                          ug=np.concatenate([g[2] for g in gas_parts]),
                          mass=np.concatenate([g[3] for g in gas_parts]))

    if nbody_components and p.ncomp > 0 and p.npart_comp[0] > 0:
        nb = nbody.build_nbody(p, verbose=False)
        for i, c in enumerate(nbody.sample_nbody(p, nb, rng, verbose=False)):
            species[i + 1] = dict(pos=c["pos"], vel=c["vel"], mass=c["mpart"])
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)
    if not species:
        raise ValueError(f"system {name!r} has no particles (check ndisc, nhalo, npart_comp)")
    return System(name, species, meta)


def read_system(path, name=None) -> System:
    """Load a GADGET HDF5 file (ours or the Fortran codes') as a System."""
    species = {}
    with h5py.File(path, "r") as f:
        masstab = np.asarray(f["Header"].attrs["MassTable"], dtype=float)
        for t in range(6):
            key = f"PartType{t}"
            if key not in f:
                continue
            g = f[key]
            s = dict(pos=g["Coordinates"][:].astype(np.float64), vel=g["Velocities"][:].astype(np.float64))
            if "Masses" in g:
                s["mass"] = g["Masses"][:].astype(np.float64)
            elif masstab[t] > 0:
                s["mass"] = float(masstab[t])
            else:
                raise ValueError(f"{path}: PartType{t} has neither Masses nor a MassTable entry")
            if "InternalEnergy" in g:
                s["ug"] = g["InternalEnergy"][:].astype(np.float64)
            species[t] = s
    return System(name or Path(path).stem, species, {"file": str(path)})
