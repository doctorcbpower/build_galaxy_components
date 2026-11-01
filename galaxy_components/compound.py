"""Combining the gas disc and the gas halo.

The Fortran workflow builds the two gas components independently and
``merge_comps.exe`` just concatenates them.  That leaves two problems where
they overlap:

  1. The hot halo gas fills the volume the disc occupies (mass is double
     counted there, and the halo particles sit inside the disc).
  2. The disc solution assumes zero external pressure, but the halo exerts a
     pressure P_h(r) = (gamma-1) rho_h u_h on it, so the two are not in
     pressure balance at the interface.

``compound_gas`` fixes both with a pressure-balance criterion.  At a point
(R, z) the cold disc has pressure P_d = rho_d cs^2 and the hot halo has
P_h(sqrt(R^2+z^2)).  The disc is taken to exist where P_d >= P_h:

  * the disc is truncated at the height where P_d falls to P_h (its vertical
    mass profile is cut and renormalised, rings with P_d(z=0) < P_h vanish);
  * halo particles are excised wherever P_d >= P_h.

What this does *not* do: the halo's hydrostatic solution ignores the disc's
gravity, and the disc's vertical solution ignores the halo's weight, so the
combination is in equilibrium only approximately (the dark-matter-only halo
potential dominates in the outskirts, where the criterion bites).  Use
``consistent=True`` in ``build_gas_halo`` so both see the same dark halo.
"""
from dataclasses import replace

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .constants import PI
from .gas_disc import GasDisc, sigma_shape
from .gas_halo import HaloGas


def _safe_pressure(halo: HaloGas, r):
    with np.errstate(invalid="ignore", divide="ignore"):
        p = halo.pressure(r)
    return np.where(np.isfinite(p), p, 0.0)


def _halo_pressure(halo, R, z):
    """Halo pressure at (R, z): the real axisymmetric one if the halo was adapted to the disc."""
    if getattr(halo, "pressure_rz", None) is not None:
        return halo.pressure_rz(R, z)
    return _safe_pressure(halo, np.hypot(R, z))


def truncate_disc(disc: GasDisc, halo: HaloGas):
    """Return a copy of ``disc`` whose vertical/radial mass tables are cut where
    the disc pressure drops to the halo pressure, plus the per-ring mass
    fraction retained."""
    R = disc.rtab[:, None]
    z = disc.ztab[None, :]
    pd = disc.cs ** 2 * 10 ** disc.lrho
    ph = _halo_pressure(halo, R, z)
    below = pd < ph
    first_below = np.where(below.any(1), below.argmax(1), disc.ztab.size)
    jc = first_below - 1                      # last z index still disc-dominated, -1 if none
    mfrac = np.where(jc >= 0, disc.mztab[np.arange(R.size), np.clip(jc, 0, None)], 0.0)
    mz = np.ones_like(disc.mztab)
    for i in range(R.size):
        if jc[i] >= 0 and mfrac[i] > 0:
            mz[i, : jc[i] + 1] = disc.mztab[i, : jc[i] + 1] / mfrac[i]
    shape = sigma_shape(disc.rtab, disc.model.rdisc_hole, disc.model.rdisc_ref, disc.model.rdisc)
    from scipy import integrate
    mr = integrate.cumulative_trapezoid(2 * PI * disc.sigma0 * disc.rtab * shape * mfrac,
                                        disc.rtab, initial=0.0)
    out = replace(disc, mztab=mz, mrtab=mr)
    return out, mfrac, jc


def disc_pressure_interpolator(disc: GasDisc):
    """P_d(R, |z|) by bilinear interpolation of log10 rho in (ln R, ln z); 0 beyond the table."""
    f = RegularGridInterpolator((np.log(disc.rtab), np.log(disc.ztab)), disc.lrho,
                                bounds_error=False, fill_value=-300.0)

    def pd(R, z):
        pts = np.column_stack((np.log(np.clip(R, disc.rtab[0], None)),
                               np.log(np.clip(np.abs(z), disc.ztab[0], None))))
        out = disc.cs ** 2 * 10 ** f(pts)
        return np.where(np.asarray(R) > disc.rtab[-1], 0.0, out)
    return pd


def excise_halo(pos, disc: GasDisc, halo: HaloGas):
    """Boolean mask of halo particles to KEEP (those not inside the disc)."""
    pos = np.asarray(pos, dtype=float)
    R = np.hypot(pos[:, 0], pos[:, 1])
    z = pos[:, 2]
    pd = disc_pressure_interpolator(disc)(R, z)
    ph = _halo_pressure(halo, R, z)
    return ~(pd >= ph)


def mass_summary(disc, trunc, mfrac):
    full = disc.mrtab[-1]
    cut = trunc.mrtab[-1]
    return {"disc_mass_full": full, "disc_mass_after_truncation": cut,
            "disc_outer_radius": float(disc.rtab[np.where(mfrac > 0)[0].max()]) if (mfrac > 0).any() else 0.0}
