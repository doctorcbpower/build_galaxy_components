"""Equilibrium N-body components via Eddington inversion (port of mk_nbody_comps.f90).

This is the default N-body path.  It is spherical and isotropic and sees only
the N-body components themselves (halo, bulge, ...): a stellar or gas disc is
NOT sampled and its gravity is NOT in the potential, so the result is exact in
isolation and only perturbatively consistent once a disc is added.

Spherical, isotropic, multi-component.  The tabulated relative potential
Psi(r) = (rscale/mtot) int_r^rmax M(r')/r'^2 dr' is built from the summed
mass profile; each component's f(E) is then obtained from

    f(E) = 1/(sqrt(8) pi^2) int_0^E (d^2 rho / dPsi^2) / sqrt(E - Psi) dPsi.

Quadratures that were Romberg/midpoint calls in Fortran are done exactly on the
piecewise-linear / log-linear interpolants the Fortran integrands implied.

Velocity anisotropy (per component, in the parameter file):

  * ``raniso_comp = r_a / rs`` -- Osipkov-Merritt, beta(r) = r^2/(r^2 + r_a^2): isotropic in
    the centre, radially biased outside.  f depends on Q = E - L^2/(2 r_a^2) and is found by
    the same Eddington inversion applied to rho_Q = (1 + r^2/r_a^2) rho.
  * ``beta_comp = beta`` -- constant beta in (-1/2, 1/2), radial or tangential:
    f(E, L) = L^(-2 beta) g(E), with g from a generalised Abel inversion of r^(2 beta) rho,

        g(E) = [K Gamma(a) Gamma(g)]^-1 int_0^E (d^2 R/dPsi^2) (E-Psi)^(g-1) dPsi,
        R = r^(2 beta) rho,  a = 3/2 - beta,  g = 1/2 + beta,  K = 2 pi 2^(1/2-beta) sqrt(pi) G(1-beta)/G(3/2-beta)

    (beta = 0 reduces to Eddington's formula).  Surface terms at Psi = 0 are dropped for the
    constant-beta model; they vanish for a tapered or r^-3 halo.
  At most one of the two may be set per component; both zero is the isotropic model.

Known Fortran quirks and what we do about them (``legacy`` flag):
  * drptab edge extrapolation mixes in dptab differences (legacy keeps this)
  * the potential at the sampled radius is extrapolated from the wrong
    segment (i, i+1) instead of (i-1, i)  (we interpolate correctly)
  * velocities are scaled with sqrt(G mass_comp(1)/rscale_comp(1)) although
    Psi is normalised by mtot (identical for NFW/Burkert; we use mtot)
  * `rs` in the f(E) prefactor is the *last* component's scale radius (we use
    component 1's, consistent with the Psi normalisation; only a constant)
"""
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.special import gamma as Gamma

from .constants import GGDT, PI, virial_mass, virial_radius
from .model import GalaxyModel
from .params import Params

NTAB = 1000
FRCUT = 1.05
PROFILES = {"nfw": 1, "hernquist": 2, "burkert": 3}


# ---------------------------------------------------------------- profiles
def rho_profile(kind, r, mass, rs, r200, galic=True):
    """Density of an (untruncated) profile with total/virial mass `mass`."""
    x = np.asarray(r, dtype=float) / rs
    rho0 = get_rho0(kind, mass, r200, rs)
    if kind == "nfw":
        return rho0 / (x * (1 + x) ** 2)
    if kind == "hernquist":
        if galic:  # rho_scale = mhalo / (2 pi a^3)  (rhodm with GalIC)
            return mass / (2 * PI * rs ** 3) / (x * (1 + x) ** 3)
        return (mass / 2 / PI / rs ** 3) * ((1 + r200 / rs) ** 2) / (r200 / rs) ** 2 / (x * (1 + x) ** 3)
    if kind == "burkert":
        return rho0 / ((1 + x) * (1 + x * x))
    raise NotImplementedError(f"profile {kind!r}")


def get_rho0(kind, mref, rref, rscale):
    x = rref / rscale
    if kind == "nfw":
        return mref / (np.log(1 + x) - x / (1 + x)) / 4 / PI / rscale ** 3
    if kind == "hernquist":
        return mref * (1 + x) ** 2 / x ** 2 / 2 / PI / rscale ** 3
    if kind == "burkert":
        return mref / (np.log(1 + x * x) + 2 * np.log(1 + x) - 2 * np.arctan(x)) / PI / rscale ** 3
    raise NotImplementedError(f"profile {kind!r}")


def mhenc(kind, r, mhalo, rs, c200, r200, galic=True):
    """Fortran mhenc (note: for Hernquist+GalIC the scale is rs*sqrt(2 f(c)))."""
    r = np.asarray(r, dtype=float)
    if kind == "nfw":
        x = r / rs
        return mhalo / (np.log(1 + c200) - c200 / (1 + c200)) * (np.log(1 + x) - x / (1 + x))
    if kind == "hernquist":
        a = rs * (np.sqrt(2 * (np.log(1 + c200) - c200 / (1 + c200))) if galic else 1.0)
        x = r / a
        m = mhalo * x * x / (1 + x) ** 2
        return m if galic else m * (1 + a / r200) ** 2
    if kind == "burkert":
        x = r / rs
        fofc = np.log(1 + c200 ** 2) + 2 * np.log(1 + c200) - 2 * np.arctan(c200)
        return mhalo / fofc * (np.log(1 + x * x) + 2 * np.log(1 + x) - 2 * np.arctan(x))
    raise NotImplementedError(f"profile {kind!r}")


def slope(kind, r, rs):
    x = np.asarray(r, dtype=float) / rs
    if kind == "nfw":
        return -(1 + 3 * x) / (1 + x)
    if kind == "hernquist":
        return -(1 + 4 * x) / (1 + x)
    return -(x / (1 + x) + 2 * x * x / (1 + x * x))


def rho_trunc(kind, r, rs, rho_scale, rt, rstar, rho_star):
    """Density tapered exponentially beyond rt (rho_trunc in structure.f90)."""
    r = np.asarray(r, dtype=float)
    x = r / rs
    if kind == "nfw":
        inner = rho_scale / (x * (1 + x) ** 2)
    elif kind == "hernquist":
        inner = rho_scale / (x * (1 + x) ** 3)
    else:
        inner = rho_scale / ((1 + x) * (1 + x * x))
    y = r / rstar
    outer = np.where(y > 50, 0.0, rho_star / (r / rt) ** 2 * np.exp(-np.minimum(y, 50)))
    return np.where(r < rt, inner, outer)


# ----------------------------------------------------------- numerics
def _stencil(f, h):
    """4th-order first derivative on interior points, linear extrapolation at edges
    exactly as in the Fortran code (rows 1,2,n-1,n)."""
    d = np.zeros_like(f)
    d[..., 2:-2] = (-f[..., 4:] + 8 * f[..., 3:-1] - 8 * f[..., 1:-3] + f[..., :-4]) / 12.0 / h
    return _edges(d, d)


def _edges(d, ref):
    """Fortran edge fill; `ref` supplies the differences (legacy bug: ref != d)."""
    d[..., 0] = d[..., 2] - 2 * (ref[..., 3] - ref[..., 2])
    d[..., 1] = d[..., 2] - (ref[..., 3] - ref[..., 2])
    d[..., -2] = d[..., -3] + (ref[..., -3] - ref[..., -4])
    d[..., -1] = d[..., -3] + 2 * (ref[..., -3] - ref[..., -4])
    return d


def _cum_logvar(y, f):
    """Cumulative int f(y) e^{3y} dy on a grid, f piecewise linear in y (exact)."""
    h = np.diff(y)
    f0, f1 = f[:-1], f[1:]
    e3 = np.exp(3 * h)
    seg = np.exp(3 * y[:-1]) * (f0 * (e3 - 1) / 3 + (f1 - f0) / h * (h * e3 / 3 - (e3 - 1) / 9))
    return np.concatenate(([0.0], np.cumsum(seg)))


def _tail_potential(rtab, mtab):
    """I_i = int_{ln r_i}^{ln r_max} M/r dln r, with ln(M/r^2) linear in ln r (as fpot)."""
    y = np.log(rtab)
    q = np.log(mtab) - 2 * y            # ln(M/r^2)
    h = np.diff(y)
    s = np.diff(q) / h                   # slope
    # segment: int exp(q0 + s(y-y0)) e^{y} dy = e^{q0+y0} (e^{(1+s)h}-1)/(1+s)
    a = 1.0 + s
    small = np.abs(a) < 1e-12
    seg = np.where(small, np.exp(q[:-1] + y[:-1]) * h,
                   np.exp(q[:-1] + y[:-1]) * np.expm1(np.where(small, 1.0, a) * h) / np.where(small, 1.0, a))
    cum_from_top = np.concatenate((np.cumsum(seg[::-1])[::-1], [0.0]))
    return cum_from_top


def _fabel_integral(etab, ptab_asc, g_asc):
    """For each E: int_{0}^{E} g(Psi)/sqrt(E-Psi) dPsi, g piecewise linear (exact)."""
    p0, p1 = ptab_asc[:-1], ptab_asc[1:]
    g0, g1 = g_asc[:-1], g_asc[1:]
    b = (g1 - g0) / (p1 - p0)
    a = g0 - b * p0
    out = np.empty(etab.size)
    for n, E in enumerate(etab):
        # segments overlapping [0,E]
        hi = np.minimum(p1, E)
        lo = np.minimum(p0, E)
        ok = hi > lo
        u_hi = E - lo[ok]
        u_lo = E - hi[ok]
        c = a[ok] + b[ok] * E
        H = lambda u: -2 * c * np.sqrt(u) + (2 * b[ok] / 3) * u ** 1.5
        out[n] = np.sum(H(u_lo) - H(u_hi))
    return out


def _abel_integral(etab, ptab_asc, g_asc, gam):
    """For each E: int_0^E g(Psi) (E-Psi)^(gam-1) dPsi with g piecewise linear (exact), 0 < gam < 1."""
    p0, p1 = ptab_asc[:-1], ptab_asc[1:]
    g0, g1 = g_asc[:-1], g_asc[1:]
    b = (g1 - g0) / (p1 - p0)
    a = g0 - b * p0
    out = np.empty(etab.size)
    for n, E in enumerate(etab):
        hi = np.minimum(p1, E)
        lo = np.minimum(p0, E)
        ok = hi > lo
        u_hi, u_lo = E - lo[ok], E - hi[ok]
        c = a[ok] + b[ok] * E
        G = lambda u: c * u ** gam / gam - b[ok] * u ** (gam + 1) / (gam + 1)
        out[n] = np.sum(G(u_hi) - G(u_lo))
    return out


def beta_prefactor(beta):
    """1 / (K Gamma(a) Gamma(gamma)) for the constant-beta inversion (1/(sqrt(8) pi^2) at beta = 0)."""
    B = np.sqrt(PI) * Gamma(1 - beta) / Gamma(1.5 - beta)
    K = 2 * PI * B * 2 ** (0.5 - beta)
    return 1.0 / (K * Gamma(1.5 - beta) * Gamma(0.5 + beta))


def clipping_effect(etab, f, expo, menc_row, ptab):
    """Largest fractional change of rho(Psi) if the negative part of f were dropped.

    rho(Psi) is rebuilt from f+ and f- separately (kernel (Psi-E)^(expo-1)); only radii holding all
    but 1e-5 of the component's mass are considered.  0 means f >= 0 where it matters.
    """
    pos = _abel_integral(etab, etab, np.where(f > 0, f, 0.0), expo)
    neg = _abel_integral(etab, etab, np.where(f < 0, -f, 0.0), expo)
    mfrac = menc_row / menc_row[-1]
    lo, hi = np.interp(1.0 - 1e-5, mfrac, ptab), np.interp(1e-5, mfrac, ptab)
    rel = (etab >= lo) & (etab <= hi) & (pos > 0)
    return float(np.max(neg[rel] / pos[rel])) if rel.any() else 0.0


@dataclass
class NBodyModel:
    rtab: np.ndarray
    mtab: np.ndarray
    menc_tab: np.ndarray    # (ncomp, ntab)
    ptab: np.ndarray        # relative potential, dimensionless
    etab: np.ndarray
    ftab: np.ndarray        # (ncomp, nE)
    mass_comp: np.ndarray
    rscale_comp: np.ndarray
    mtot_norm: float
    vscale: float           # sqrt(G mtot / rscale1)  [km/s]
    r200: float
    m200: float
    kinds: tuple = ()       # per component: "iso", "om" or "beta"
    aniso: tuple = ()       # per component: r_a [kpc] (om), beta (beta), 0 (iso)


def build_nbody(p: Params, verbose=True, legacy=True) -> NBodyModel:
    log = print if verbose else (lambda *a, **k: None)
    ncomp = p.ncomp
    if p.istaper and any(p.rtrunc_comp[i] <= 0 for i in range(ncomp)):
        raise ValueError("istaper needs rtrunc_comp > 0 for every component "
                         f"(got {list(p.rtrunc_comp[:ncomp])})")
    if p.mdisc != 0 or p.mstar != 0:
        log("Note: disc mass in &DISC is not included in the spherical N-body potential; "
            "the disc is gas-only and acts only as a perturbation on these components.")
    m200 = virial_mass(p.v200) if p.v200 > 0 else p.m200
    r200 = virial_radius(m200)
    rs = r200 / p.c200
    log(f"Virial radius [kpc/h]: {r200:.6g}\nScale radius [kpc/h]: {rs:.6g}")

    kinds = [p.form_comp[i] or "nfw" for i in range(ncomp)]
    for k in kinds:
        if k not in PROFILES:
            raise NotImplementedError(f"form_comp {k!r}: only nfw/hernquist/burkert are ported")
    mass_comp = np.array([p.fmass_comp[i] * m200 for i in range(ncomp)])
    rscale = np.array([p.rscale_comp[i] * rs for i in range(ncomp)])
    rtrunc = np.array([p.rtrunc_comp[i] * r200 for i in range(ncomp)])

    lrmin, lrmax = -4.0 * np.log(10), 4.0 * np.log(10)
    ly = np.linspace(lrmin, lrmax, NTAB)
    dlr = ly[1] - ly[0]
    rtab = np.exp(ly)

    # taper parameters per component (taper_mass_profile)
    rstar = np.zeros(ncomp)
    rhostar = np.zeros(ncomp)
    rho_comp = np.zeros(ncomp)
    rhotab = np.zeros((ncomp, NTAB))
    for i, k in enumerate(kinds):
        rho_comp[i] = get_rho0(k, mass_comp[i], r200, rscale[i])
        if p.istaper:
            # The Fortran taper_mass_profile runs before rscale_comp is applied and so uses the global
            # NFW scale rs = r200/c200 for the slope and density at the truncation radius.  For a
            # component with a different scale that makes its density discontinuous at rtrunc and its
            # distribution function negative, so use the component's own scale.  (Identical for a single
            # component with rscale_comp = 1, which is what the Fortran comparison covers.)
            sl = slope(k, rtrunc[i], rscale[i])
            rstar[i] = -rtrunc[i] / (2.0 + sl)
            rho_at_rt = rho_profile(k, rtrunc[i], mass_comp[i], rscale[i], r200, galic=True)
            rhostar[i] = rho_at_rt * np.exp(-2.0 - sl)
            rhotab[i] = rho_trunc(k, rtab, rscale[i], rho_comp[i], rtrunc[i], rstar[i], rhostar[i])
        else:
            rhotab[i] = rho_profile(k, rtab, mass_comp[i], rscale[i], r200, galic=True)

    # enclosed masses
    menc_tab = np.zeros((ncomp, NTAB))
    mtab = np.zeros(NTAB)
    for i, k in enumerate(kinds):
        integrand = 4 * PI * rhotab[i]
        cum = _cum_logvar(ly, integrand)
        if p.istaper:  # mass frozen beyond frcut * rtrunc
            cut = rtab >= FRCUT * rtrunc[i]
            if cut.any():
                cum[cut] = cum[np.argmax(cut)]
        m1 = mhenc(k, rtab[0], mass_comp[i], rscale[i], p.c200, r200, galic=True)
        menc_tab[i] = m1 + cum
        mtab += menc_tab[i]
    if p.mbh > 0:
        mtab[0] += p.mbh  # (Fortran adds the BH at the first bin only)

    mtot_r200 = np.interp(r200, rtab, mtab)
    rscale1 = rscale[0]
    log(f"Total mass [1e10 Msol/h]: {mtot_r200:.8g}")
    log(f"Circular speed at virial radius [km/s]: {np.sqrt(GGDT * m200 / r200):.6g}")

    mtot = mass_comp[0] if kinds[0] in ("nfw", "burkert") else mtab[-1]

    # relative potential Psi(r) = (rscale/mtot) int_r^rmax M/r'^2 dr'
    ptab = _tail_potential(rtab, mtab) * rscale1 / mtot

    # anisotropy: the inversion is applied to an effective density (see module docstring)
    kinds, aniso = [], []
    rho_eff = rhotab.copy()
    for i in range(ncomp):
        ra, beta = p.raniso_comp[i], p.beta_comp[i]
        if ra < 0:
            raise ValueError(f"raniso_comp[{i}] must be >= 0")
        if ra > 0 and beta != 0:
            raise ValueError(f"component {i + 1}: set raniso_comp (Osipkov-Merritt) or beta_comp, not both")
        if ra > 0:
            kinds.append("om")
            aniso.append(ra * rs)
            rho_eff[i] = rhotab[i] * (1.0 + (rtab / (ra * rs)) ** 2)
        elif beta != 0:
            if not -0.5 < beta < 0.5:
                raise ValueError(f"beta_comp[{i}] = {beta}: the constant-beta inversion is implemented for "
                                 "-1/2 < beta < 1/2 (a cusp rho ~ r^-1 cannot be more radial than beta = 1/2)")
            kinds.append("beta")
            aniso.append(float(beta))
            rho_eff[i] = rhotab[i] * rtab ** (2 * beta)
        else:
            kinds.append("iso")
            aniso.append(0.0)

    # derivatives of rho wrt ln r, 1st and 2nd (4th order)
    drtab = _stencil(rho_eff, dlr)
    ddrtab = _stencil(drtab, dlr)

    dptab = rtab ** 2 / GGDT / mtab
    drp = np.zeros(NTAB)
    lm = np.log(mtab)
    drp[2:-2] = (-lm[4:] + 8 * lm[3:-1] - 8 * lm[1:-3] + lm[:-4]) / 12.0 / dlr
    drp = _edges(drp, dptab if legacy else drp)
    drp = (2.0 - drp) * rtab / GGDT / mtab

    pref = GGDT ** 2 * mtot / rscale1 ** 2
    drhopsi = ((ddrtab - drtab) / rtab ** 2 * dptab ** 2 + drtab / rtab * dptab * drp) * pref

    # distribution function on a log grid in E
    lemin, lemax = np.log10(ptab[-2]), np.log10(ptab[0])
    etab = 10 ** np.linspace(lemin, lemax, NTAB)
    # Psi decreases with r: reverse so it increases
    pasc = ptab[::-1]
    ftab = np.empty((ncomp, NTAB))
    for i in range(ncomp):
        if kinds[i] == "iso":
            raw_iso = _fabel_integral(etab, pasc, drhopsi[i][::-1]) / np.sqrt(8.0) / PI ** 2
            effect = clipping_effect(etab, raw_iso, 1.5, menc_tab[i], ptab)
            if effect > 0.01:
                warnings.warn(f"component {i + 1}: the isotropic distribution function is negative where it matters "
                              f"(clipping it to zero changes the component's density by up to {100 * effect:.2g} %), "
                              "so the sampled component will not reproduce its density profile. This usually means "
                              "a centrally concentrated or sharply truncated component sitting in a much more "
                              "massive halo.", RuntimeWarning, stacklevel=2)
            ftab[i] = raw_iso
            ftab[i][ftab[i] <= 1e-12] = 0.0
            continue
        if kinds[i] == "om":
            raw = _fabel_integral(etab, pasc, drhopsi[i][::-1]) / np.sqrt(8.0) / PI ** 2
            # surface term (1/sqrt(Q)) d rho_Q/dPsi |_0 : nonzero for an untapered r^-3 halo (rho_Q ~ Psi)
            b_code = -drtab[i, -3] * rtab[-3] / (rscale1 * mtab[-3])
            raw = raw + b_code / np.sqrt(8.0) / PI ** 2 / np.sqrt(etab)
        else:
            beta = aniso[i]
            raw = beta_prefactor(beta) * _abel_integral(etab, pasc, drhopsi[i][::-1], 0.5 + beta)
        # A negative f means the anisotropy is not realisable.  Tiny negative dips (numerical noise in the
        # second derivative, at energies that carry no mass) are harmless, so judge by the effect on the
        # density: refuse only if dropping f- would change rho by > 0.1 %.
        expo = 1.5 if kinds[i] == "om" else 1.5 - aniso[i]
        worst = clipping_effect(etab, raw, expo, menc_tab[i], ptab)
        if worst > 1e-3:
            raise ValueError(f"component {i + 1}: the {kinds[i]} distribution function is negative "
                             f"(clipping it would change the density by {100 * worst:.2g} %), so this anisotropy "
                             "is not realisable for this density profile; "
                             + ("increase raniso_comp" if kinds[i] == "om" else "reduce |beta_comp|"))
        raw[raw < 0] = 0.0
        ftab[i] = raw
        log(f"component {i + 1}: {'Osipkov-Merritt r_a = %.4g kpc' % aniso[i] if kinds[i] == 'om' else 'constant beta = %g' % aniso[i]}")

    return NBodyModel(rtab, mtab, menc_tab, ptab, etab, ftab, mass_comp, rscale, mtot,
                      np.sqrt(GGDT * mtot / rscale1), r200, m200, tuple(kinds), tuple(aniso))


def _local_frame(pos):
    """Orthonormal (e_r, e_1, e_2) at each position."""
    er = pos / np.linalg.norm(pos, axis=1)[:, None]
    ref = np.where(np.abs(er[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    e1 = ref - np.sum(ref * er, axis=1)[:, None] * er
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    return er, e1, np.cross(er, e1)


def _sample_s(table_e, table_f, phi, expo, rng, chunk=20000, sgrid=192, smin=1e-9):
    """Draw s = Psi - E with density proportional to f(Psi - s) s^expo on (smin*Psi, Psi].

    Proposals are log-uniform in s and are accepted against a per-particle maximum found on a
    log grid, so the cost stays O(1) draws however steeply f rises towards the deepest energies
    (a cusp makes f span ~12 decades, which defeats a uniform velocity-box proposal).
    """
    n = len(phi)
    out = np.empty(n)
    lg = np.log(smin)
    grid = np.exp(lg * (1.0 - (np.arange(sgrid) + 0.5) / sgrid))          # fractions of Psi, log-spaced
    f = lambda e: np.interp(e, table_e, table_f, left=0.0)
    for lo in range(0, n, chunk):
        ph = phi[lo:lo + chunk]
        sg = ph[:, None] * grid[None, :]
        hmax = 1.3 * (f(ph[:, None] - sg) * sg ** (expo + 1.0)).max(axis=1)   # x s: log-uniform proposal ~ 1/s
        if np.any(hmax <= 0):
            raise RuntimeError("distribution function vanishes at a sampled radius")
        res = np.empty(len(ph))
        todo = np.arange(len(ph))
        for _ in range(2000):
            sv = ph[todo] * np.exp(lg * rng.random(todo.size))
            h = f(ph[todo] - sv) * sv ** (expo + 1.0)
            ok = h > hmax[todo] * rng.random(todo.size)
            res[todo[ok]] = sv[ok]
            todo = todo[~ok]
            if not todo.size:
                break
        else:
            raise RuntimeError("rejection sampling did not converge (check the distribution function)")
        out[lo:lo + chunk] = res
    return out


def _velocities_om(nb, i, radius, phi, ra, rng):
    """Osipkov-Merritt velocities.  With u = v_r^2/2, w = k v_t^2/2, k = 1 + r^2/r_a^2, Q = Psi - u - w
    and d^3v = (2 pi/k) dw du / sqrt(2u):  p(Q) ~ f(Q) sqrt(Psi-Q), and for fixed s = Psi-Q the split
    u : w has density ~ u^(-1/2), i.e. u = s U^2."""
    n = len(radius)
    k = 1.0 + (radius / ra) ** 2
    s = _sample_s(nb.etab, nb.ftab[i], phi, 0.5, rng)
    u = s * rng.random(n) ** 2
    vr = np.sqrt(2.0 * u) * np.where(rng.random(n) < 0.5, -1.0, 1.0)
    vt = np.sqrt(2.0 * (s - u) / k)
    psi = 2 * PI * rng.random(n)
    return vr, vt * np.cos(psi), vt * np.sin(psi)


def _velocities_beta(nb, i, phi, beta, rng):
    """Constant beta: f d^3v ~ g(E) v^(2-2b) sin^(1-2b)(eta) dv d(eta), so the angle eta from the
    radial direction separates from the speed; with s = v^2/2 the speed density is g(Psi-s) s^(1/2-b)."""
    n = len(phi)
    eta_grid = np.linspace(0.0, PI, 4001)
    w = np.sin(eta_grid) ** (1.0 - 2.0 * beta)
    cdf = np.concatenate(([0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(eta_grid))))
    eta = np.interp(rng.random(n), cdf / cdf[-1], eta_grid)
    speed = np.sqrt(2.0 * _sample_s(nb.etab, nb.ftab[i], phi, 0.5 - beta, rng))
    return speed * np.cos(eta), speed * np.sin(eta)


def sample_nbody(p: Params, nb: NBodyModel, rng, verbose=True):
    """Sample particles.  Returns per-component dict list: pos, vel (float32), mpart."""
    log = print if verbose else (lambda *a, **k: None)
    out = []
    ptab_asc, rtab_by_p = nb.ptab[::-1], nb.rtab[::-1]
    for i in range(p.ncomp):
        mpart = nb.mass_comp[i] / p.npart_comp[i]
        mtot_c = nb.menc_tab[i, -1]
        n = 1 + int(np.floor(mtot_c / mpart))
        # radii: invert the cumulative mass
        radius = np.interp(rng.random(n) * mtot_c, nb.menc_tab[i], nb.rtab)
        phi = np.interp(radius, nb.rtab, nb.ptab)
        kind = nb.kinds[i] if nb.kinds else "iso"
        vel = np.zeros((n, 3))
        if kind == "iso":
            fmax = np.interp(phi, nb.etab, nb.ftab[i])  # f at E=Psi (v=0)
            todo = np.arange(n)
            attempts = np.zeros(n, dtype=int)
            while todo.size:
                vesc = np.sqrt(2.0 * phi[todo])
                v = (2.0 * rng.random((todo.size, 3)) - 1.0) * vesc[:, None]
                e = phi[todo] - 0.5 * np.sum(v * v, axis=1)
                ok = (e > 0) & (e >= nb.etab[0]) & (e <= nb.etab[-1])
                fs = np.interp(e, nb.etab, nb.ftab[i])
                ok &= fs > fmax[todo] * rng.random(todo.size)
                vel[todo[ok]] = v[ok]
                attempts[todo] += 1
                todo = todo[~ok]
                # Fortran restarts with a new radius after 1000 failed attempts
                redo = todo[attempts[todo] > 1000]
                if redo.size:
                    radius[redo] = np.interp(rng.random(redo.size) * mtot_c, nb.menc_tab[i], nb.rtab)
                    phi[redo] = np.interp(radius[redo], nb.rtab, nb.ptab)
                    fmax[redo] = np.interp(phi[redo], nb.etab, nb.ftab[i])
                    attempts[redo] = 0
        if kind != "iso":
            for _ in range(50):
                bad = np.interp(phi, nb.etab, nb.ftab[i]) <= 0
                if not bad.any():
                    break
                radius[bad] = np.interp(rng.random(bad.sum()) * mtot_c, nb.menc_tab[i], nb.rtab)
                phi[bad] = np.interp(radius[bad], nb.rtab, nb.ptab)
        # isotropic directions
        x1 = 2.0 * rng.random(n) - 1.0
        az = 2 * PI * rng.random(n)
        rc = np.sqrt(1.0 - x1 * x1)
        pos = radius[:, None] * np.column_stack((x1, rc * np.cos(az), rc * np.sin(az)))
        if kind != "iso":      # velocities in the local frame (e_r, e_1, e_2) at each position
            er, e1, e2 = _local_frame(pos)
            if kind == "om":
                vr, v1, v2 = _velocities_om(nb, i, radius, phi, nb.aniso[i], rng)
            else:
                vr, vt = _velocities_beta(nb, i, phi, nb.aniso[i], rng)
                psi = 2 * PI * rng.random(n)
                v1, v2 = vt * np.cos(psi), vt * np.sin(psi)
            vel = vr[:, None] * er + v1[:, None] * e1 + v2[:, None] * e2
        out.append({"pos": pos.astype(np.float32),
                    "vel": (vel * nb.vscale).astype(np.float32),
                    "mpart": mpart, "radius": radius})
        log(f"component {i + 1}: {n} particles of mass {mpart:.4e}")
    return out


def summarise(nb, comps):
    pos = np.concatenate([c["pos"] for c in comps]).astype(float)
    vel = np.concatenate([c["vel"] for c in comps]).astype(float)
    r = np.linalg.norm(pos, axis=1)
    sel = r <= nb.r200
    v, x = vel[sel], pos[sel]
    vr = np.sum(v * x, axis=1) / r[sel]
    return {"nvir": int(sel.sum()), "com": x.mean(0), "vcom": v.mean(0),
            "sigma": v.std(0), "sigma_tot": np.sqrt(v.var(0).sum()),
            "vr": vr.mean(), "sigma_r": vr.std()}


def measure_beta(pos, vel, edges):
    """Anisotropy beta = 1 - (s_theta^2 + s_phi^2)/(2 s_r^2) in radial bins.

    Returns (r_mid, beta, n_per_bin); second moments about zero (the halo has no mean flow).
    """
    pos = np.asarray(pos, dtype=float)
    vel = np.asarray(vel, dtype=float)
    r = np.linalg.norm(pos, axis=1)
    vr = np.sum(vel * pos, axis=1) / r
    vt2 = np.sum(vel ** 2, axis=1) - vr ** 2
    idx = np.digitize(r, edges) - 1
    mid, beta, cnt = [], [], []
    for b in range(len(edges) - 1):
        sel = idx == b
        if sel.sum() < 20:
            mid.append(np.nan), beta.append(np.nan), cnt.append(int(sel.sum()))
            continue
        mid.append(np.exp(np.mean(np.log(r[sel])))), cnt.append(int(sel.sum()))
        beta.append(1.0 - 0.5 * np.mean(vt2[sel]) / np.mean(vr[sel] ** 2))
    return np.array(mid), np.array(beta), np.array(cnt)
