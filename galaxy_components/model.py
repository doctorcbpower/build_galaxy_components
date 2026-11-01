"""Derive the halo/bulge/disc parameters from the namelist and build components."""
from dataclasses import dataclass, field

import numpy as np

from .constants import (GGDT, RHOCRIT0, DELTAVIR, virial_mass, virial_radius)
from .params import Params
from .potentials import (Component, ExponentialDisc, Hernquist, LogarithmicHalo, NFW)


@dataclass
class GalaxyModel:
    m200: float
    r200: float
    rs: float
    c200: float
    mhalo: float
    mbulge: float
    mdisc: float
    mstar: float
    mbh: float
    rdisc: float
    rbulge: float
    rdisc_hole: float
    rdisc_ref: float
    halo_a: float  # Hernquist scale length (== rs for NFW)
    components: list = field(default_factory=list)  # fixed background (halo, bulge, stellar disc)

    def vc2_background(self, R):
        return sum(c.vc2(R) for c in self.components)

    def phi_background(self, R, z):
        return sum(c.phi(R[:, None], z[None, :]) if not isinstance(c, ExponentialDisc)
                   else c.phi(R, z) for c in self.components)


def build_model(p: Params, verbose=True) -> GalaxyModel:
    log = print if verbose else (lambda *a, **k: None)
    m200 = p.m200
    if p.v200 > 0.0:
        m200 = virial_mass(p.v200)
    r200 = virial_radius(m200)
    rs = r200 / p.c200
    fofc = np.log(1 + p.c200) - p.c200 / (1 + p.c200)
    hq_scale = np.sqrt(2 * fofc)  # Hernquist a = rs * sqrt(2 f(c)) (matches NFW inner profile)

    # Masses: values > 0 are fractions of M200; values <= 0 are absolute (1e10 Msol)
    mbulge = p.mbulge * m200 if p.mbulge > 0 else abs(p.mbulge)
    mdisc = p.mdisc * m200 if p.mdisc > 0 else abs(p.mdisc)
    mstar = p.mstar * mdisc if p.mstar > 0 else abs(p.mstar)
    mbh = p.mbh
    mhalo = m200 - mbulge - mdisc - mbh

    if p.halo_type == "hernquist":
        a = rs * hq_scale
        rdisc, rbulge = p.fdisc * a, p.fbulge * a
        rhole = p.fhole * rdisc
        rref = 0.0
    elif p.halo_type == "nfw":
        a = rs
        rdisc, rbulge = p.fdisc * rs, p.fbulge * rs
        rhole = p.fhole * rdisc
        rref = 0.0
    elif p.halo_type == "logarithmic":
        a = rs
        rdisc, rbulge = abs(p.fdisc), abs(p.fbulge)
        rhole, rref = abs(p.fhole), abs(p.fdref)
    else:
        raise ValueError(f"unknown halo_type {p.halo_type!r}")

    comps: list[Component] = []
    if p.halo_type == "hernquist":
        mtot = mhalo if p.galic else mhalo * (1 + a / r200) ** 2
        comps.append(Hernquist(mtot, a))
    elif p.halo_type == "nfw":
        comps.append(NFW(mhalo, rs, p.c200))
    else:
        comps.append(LogarithmicHalo(p.v200, p.rcore, p.qflat))
    if mbulge > 0:
        comps.append(Hernquist(mbulge, rbulge))
    if mstar > 0:
        comps.append(ExponentialDisc(mstar, rdisc, z0=p.zdisc * rdisc))

    log(f"Halo ({p.halo_type}): M200={m200:g} R200={r200:g} rs={rs:g} a={a:g}")
    return GalaxyModel(m200, r200, rs, p.c200, mhalo, mbulge, mdisc, mstar, mbh,
                       rdisc, rbulge, rhole, rref, a, comps)
