"""Hot gaseous halo in hydrostatic equilibrium (port of mk_gas_halo.f90 / gas_halo.f90).

Polytropic gas, P = A rho^Gamma, in the potential of a dark halo, following
Komatsu & Seljak (2001).  With f(x) the normalised enclosed mass of the halo
(x = r/a),

    rho_g(x) = rho0 [1 - (3/eta0) (1-1/Gamma) (c/f(c)) int_0^x f(u)/u^2 du]^(1/(Gamma-1))
    u_g(x)   = u0 (rho_g/rho0)^(Gamma-1)

where Gamma is set by matching the gas and dark slopes at r200 (K&S eq. 25) and
eta0 by the same matching.  Both integrals have closed forms for NFW and
Hernquist, so no quadrature is needed except for the gas mass profile.

``consistent=True`` (default) makes the halo potential the *same* dark halo the
disc and N-body codes use (Hernquist scale length a = rs*sqrt(2 f(c)) under
GalIC, and M_halo(<r200) rather than M200 as the mass normalisation).
``consistent=False`` reproduces the Fortran, which uses a = r200/c200 and M200
regardless -- a different potential from the one the disc sees.
"""
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import integrate

from .constants import GGDT, KB, KMS_IN_MS, MP, PI, GAMMA, virial_mass, virial_radius
from .params import Params

MU_HALO = 0.61  # mean molecular weight of the ionised halo gas
U_TO_TEMP = ((GAMMA - 1) * MP / KB) * KMS_IN_MS ** 2


# ---- dark-halo helper functions (normalised: x = r/a) ----------------------
def _fnorm(kind, x):
    x = np.asarray(x, dtype=float)
    if kind == "nfw":
        return np.log1p(x) - x / (1 + x)
    return x * x / (1 + x) ** 2


def _int_f_over_r2(kind, a, b):
    """int_a^b f(u)/u^2 du in closed form."""
    if kind == "nfw":
        g = lambda t: -np.log1p(t) / t
    else:
        g = lambda t: -1.0 / (1 + t)
    return g(b) - g(a)


def _slope(kind, x):
    return -(1 + 3 * x) / (1 + x) if kind == "nfw" else -(1 + 4 * x) / (1 + x)


@dataclass
class HaloGas:
    kind: str
    a: float                 # dark-halo scale length used [kpc]
    m200: float
    r200: float
    conc: float              # r200 / a
    gamma_poly: float
    a0: float
    temp0: float             # central temperature [K]
    rho0: float              # central gas density [1e10 Msol/kpc^3]
    mgas: float              # gas mass inside r200
    rtab: np.ndarray         # [kpc]
    mtab: np.ndarray         # M_gas(<r)/mgas
    utab: np.ndarray         # specific internal energy [(km/s)^2]
    rho: Callable            # rho_gas(r)
    u: Callable              # u_gas(r)
    m_ref: float = 0.0       # dark-halo mass inside r200 used in the potential
    mdm: Callable = None     # enclosed gravitating mass M(<r) [1e10 Msol] the gas is in equilibrium with
    redge: float = 0.0       # radius where P = 0 [kpc]; 0 for a free-edge halo

    def pressure(self, r):
        return (GAMMA - 1) * self.rho(r) * self.u(r)


def build_gas_halo(p: Params, ntab=1000, consistent=True, verbose=True) -> HaloGas:
    log = print if verbose else (lambda *a, **k: None)
    kind = "hernquist" if (p.form or "nfw").lower() == "hernquist" else "nfw"
    m200 = virial_mass(p.v200) if p.v200 > 0 else p.m200
    r200 = virial_radius(m200)
    rs = r200 / p.c200
    fofc = np.log(1 + p.c200) - p.c200 / (1 + p.c200)
    if consistent and kind == "hernquist" and p.galic:
        a = rs * np.sqrt(2 * fofc)
        conc = r200 / a
        m_ref = m200 * float(_fnorm(kind, conc))   # M_halo(<r200) with total mass m200
    else:
        a, conc, m_ref = rs, p.c200, m200
    router = p.rtrunc * r200
    log(f"Gas halo ({kind}): M200={m200:g} R200={r200:g} a={a:g} c_eff={conc:g}")

    gam = 1.15 + 0.01 * (conc - 6.5)  # K&S eq. 25
    xmin = 1e-8                        # in units of r200
    x1, x2 = np.log(conc * xmin), np.log(conc)
    fc = float(_fnorm(kind, conc))
    a0 = 3 * (gam - 1) * (conc / fc) * _int_f_over_r2(kind, np.exp(x1), conc)
    a0 += (-3.0 / _slope(kind, conc)) * (fc / conc) * (conc / fc)
    a0 /= gam
    temp0 = (1.0 / 3.0) * (MU_HALO * MP / KB) * KMS_IN_MS ** 2 * (GGDT * m_ref / r200) * a0
    log(f"Using a gamma of {gam:.6f}; central halo temperature [K]: {temp0:.6g}")

    fac = (1.0 - 1.0 / gam) * (conc / fc)

    def rhog(x):  # lower limit 1e-6 as in gas_halo.f90
        x = np.asarray(x, dtype=float)
        return (1.0 - (3.0 / a0) * fac * _int_f_over_r2(kind, 1e-6, x)) ** (1.0 / (gam - 1.0))

    xmax = 1.1 * router / r200
    xfac = (xmax / 1.01 / xmin) ** (1.0 / (ntab - 1))
    xs = 1.01 * xmin * conc * xfac ** np.arange(ntab)       # in units of a
    # cumulative int x^2 rhog dx from exp(x1): refine the table grid for Simpson accuracy
    sub = 8
    lxs = np.log(xs)
    fine = np.unique(np.concatenate(([x1], np.interp(np.arange((ntab - 1) * sub + 1) / sub,
                                                     np.arange(ntab), lxs))))
    xf = np.exp(fine)
    cum = integrate.cumulative_simpson(xf ** 3 * rhog(xf), x=fine, initial=0.0)  # dx = x dlnx
    total = np.interp(np.log(conc), fine, cum)
    rho0 = p.fbaryon * m200 / (4 * PI * a ** 3) / total
    mgas = p.fbaryon * m200
    mtab = np.interp(lxs, fine, cum) / total
    rtab = xs * a
    utab = temp0 * rhog(xs) ** (gam - 1) / U_TO_TEMP / MU_HALO

    rho_fn = lambda r: rho0 * rhog(np.asarray(r) / a)
    u_fn = lambda r: temp0 * rhog(np.asarray(r) / a) ** (gam - 1) / U_TO_TEMP / MU_HALO
    mdm = lambda r: m_ref * _fnorm(kind, np.asarray(r, dtype=float) / a) / fc
    return HaloGas(kind, a, m200, r200, conc, gam, a0, temp0, rho0, mgas, rtab, mtab, utab, rho_fn, u_fn, m_ref,
                   mdm)


@dataclass
class DiscHaloGas(HaloGas):
    """Axisymmetric halo gas in the potential of the dark halo *plus* a fixed exponential disc.

    The gas follows the equation of state of the spherical halo it was built from, P(h), where
    h = int dP/rho is the enthalpy.  Hydrostatic equilibrium in any potential is h + Phi = const, so
    with Phi = Phi_dm(r) + Phi_disc(R, z) the enthalpy is h(R, z) = h_dm(r) - [Phi_disc(R, z) -
    Phi_disc(pole)], h_dm(r) = int_r^R_edge g_dm dr, and the density is rho(h).  The zero-pressure
    surface is the equipotential through the pole of the old edge.  The disc pulls the gas in towards
    itself; the gas mass inside r200 is renormalised (an entropy rescale, rho -> s rho, P -> s P).
    ``rho``/``u``/``pressure`` of this class are averages over angle; use ``*_rz`` for the real thing.
    """
    rho_rz: Callable = None
    pressure_rz: Callable = None
    mu_grid: np.ndarray = None      # cos(theta) grid on [0, 1]
    mu_cdf: np.ndarray = None       # (nr, nmu) conditional CDF of cos(theta) at each rtab radius
    disc: object = None

    def u_rz(self, R, z):
        rho = self.rho_rz(R, z)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(rho > 0, self.pressure_rz(R, z) / ((GAMMA - 1) * np.where(rho > 0, rho, 1.0)), 0.0)

    def draw_cos_theta(self, r, rng):
        """Signed cos(theta) for particles at radii r, from the conditional density at that radius."""
        i = np.clip(np.searchsorted(self.rtab, r) , 1, len(self.rtab) - 1)
        i = np.where(np.abs(np.log(r / self.rtab[i - 1])) < np.abs(np.log(r / self.rtab[i])), i - 1, i)
        rows = self.mu_cdf[i]
        u = rng.random(len(r))
        j = np.clip((rows < u[:, None]).sum(1), 1, rows.shape[1] - 1)
        lo, hi = rows[np.arange(len(r)), j - 1], rows[np.arange(len(r)), j]
        t = np.where(hi > lo, (u - lo) / np.where(hi > lo, hi - lo, 1.0), 0.0)
        mu = self.mu_grid[j - 1] + t * (self.mu_grid[j] - self.mu_grid[j - 1])
        return np.where(rng.random(len(r)) < 0.5, -mu, mu)


def adapt_to_disc(halo: HaloGas, disc, nmu=240, ngrid=200, verbose=True) -> DiscHaloGas:
    """Put the gas of a zero-pressure-edge halo into the potential of halo + ``disc`` (an
    ``ExponentialDisc``); see ``DiscHaloGas``."""
    from scipy.interpolate import RegularGridInterpolator
    log = print if verbose else (lambda *a, **k: None)
    if halo.redge <= 0 or halo.mdm is None:
        raise ValueError("halo_feels_disc needs a gas halo with a zero-pressure edge: set halo_zero_edge = true "
                         "with rtrunc > 1, or use a truncated N-body halo (istaper)")
    r, lr, redge = halo.rtab, np.log(halo.rtab), halo.redge
    g = GGDT * halo.mdm(r) / r ** 2
    cum = integrate.cumulative_simpson(g * r, x=lr, initial=0.0)
    hdm = np.maximum(np.interp(np.log(redge), lr, cum) - cum, 0.0)          # int_r^redge g dr
    hdm_of_r = lambda rr: np.interp(np.log(np.maximum(rr, r[0])), lr, hdm, right=0.0)
    # equation of state of the spherical solution, ascending in h
    keep = (r > 1e-5 * halo.r200) & (r < redge)
    h_tab = hdm[keep][::-1]
    rho_tab = halo.rho(r[keep])[::-1]
    p_tab = ((GAMMA - 1) * halo.rho(r[keep]) * halo.u(r[keep]))[::-1]
    h_tab, idx = np.unique(h_tab, return_index=True)
    rho_tab, p_tab = rho_tab[idx], p_tab[idx]
    h_tab, rho_tab, p_tab = np.concatenate(([0.0], h_tab)), np.concatenate(([0.0], rho_tab)), np.concatenate(([0.0], p_tab))
    hc = h_tab[-1]

    def eos(h, tab):
        h = np.asarray(h, dtype=float)
        top = tab[-1] * (np.maximum(h, 1e-300) / hc) ** (1.5 if tab is rho_tab else 2.5)   # constant-K continuation
        return np.where(h <= 0, 0.0, np.where(h <= hc, np.interp(h, h_tab, tab), top))

    # disc potential relative to the pole of the old edge, on a log grid in (R, |z|)
    lo, hi = 0.02 * disc.rd, 1.3 * redge
    Rg = np.geomspace(lo, hi, ngrid)
    phid = disc.phi_grid(Rg, Rg)
    pole = float(disc.phi_multipole(1e-9 * redge, redge))
    interp = RegularGridInterpolator((np.log(Rg), np.log(Rg)), phid - pole, bounds_error=False, fill_value=None)

    def dphi(R, z):
        R, z = np.broadcast_arrays(np.asarray(R, dtype=float), np.asarray(z, dtype=float))
        pts = np.stack((np.log(np.clip(R, lo, hi)), np.log(np.clip(np.abs(z), lo, hi))), -1)
        return interp(pts)

    def h_of(R, z):
        return hdm_of_r(np.hypot(R, z)) - dphi(R, z)

    # mass normalisation inside r200, and the grids used for sampling
    nr = 600
    rg = np.geomspace(1e-4 * halo.r200, redge, nr)
    mu = np.linspace(0.0, 1.0, nmu)
    RR = rg[:, None] * np.sqrt(1 - mu ** 2)[None, :]
    ZZ = rg[:, None] * mu[None, :]
    rho0 = eos(h_of(RR, ZZ), rho_tab)
    shell = 4 * PI * rg ** 3 * np.trapezoid(rho0, mu, axis=1)               # dM/dlnr
    menc = integrate.cumulative_simpson(shell, x=np.log(rg), initial=0.0)
    s = halo.mgas / float(np.interp(halo.r200, rg, menc))
    rho_tab, p_tab, menc = rho_tab * s, p_tab * s, menc * s
    rho1 = rho0 * s
    cdf = integrate.cumulative_trapezoid(rho1, mu, axis=1, initial=0.0)
    cdf = cdf / np.where(cdf[:, -1:] > 0, cdf[:, -1:], 1.0)
    rho_rz = lambda R, z: eos(h_of(R, z), rho_tab)
    p_rz = lambda R, z: eos(h_of(R, z), p_tab)
    lg = np.log(rg)
    rho_avg = np.trapezoid(rho1, mu, axis=1)
    p_avg = np.trapezoid(eos(h_of(RR, ZZ), p_tab), mu, axis=1)
    u_avg = np.where(rho_avg > 0, p_avg / ((GAMMA - 1) * np.maximum(rho_avg, 1e-300)), 0.0)
    log(f"Gas halo adapted to the disc (M={disc.mass:g}, Rd={disc.rd:g}): gas mass inside r200 renormalised by "
        f"{s:.4g}; midplane density at 5 kpc is {float(rho_rz(np.array([5.0]), np.array([0.0]))[0]) / float(halo.rho(5.0)):.3g}"
        f" x the spherical halo's")
    return DiscHaloGas(halo.kind, halo.a, halo.m200, halo.r200, halo.conc, halo.gamma_poly, halo.a0, halo.temp0,
                       halo.rho0, halo.mgas, rg, menc / halo.mgas, u_avg,
                       lambda rr: np.interp(np.log(np.asarray(rr, dtype=float)), lg, rho_avg),
                       lambda rr: np.interp(np.log(np.asarray(rr, dtype=float)), lg, u_avg),
                       halo.m_ref, halo.mdm, halo.redge, rho_rz, p_rz, mu, cdf, disc)


def _loglin_interp(rr, r, y, floor=1e-300):
    """Interpolate a positive profile log-log, but linearly in the one cell that touches a zero
    (so a profile that reaches zero at an edge does not get a fake exponential plunge there)."""
    rr = np.asarray(rr, dtype=float)
    i = np.clip(np.searchsorted(r, rr) - 1, 0, len(r) - 2)
    y0, y1 = y[i], y[i + 1]
    t = (rr - r[i]) / (r[i + 1] - r[i])
    lin = y0 + (y1 - y0) * t
    ok = (y0 > floor) & (y1 > floor)
    tl = (np.log(rr) - np.log(r[i])) / (np.log(r[i + 1]) - np.log(r[i]))
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.exp(np.log(np.where(ok, y0, 1.0)) * (1 - tl) + np.log(np.where(ok, y1, 1.0)) * tl)
    out = np.where(ok, lg, lin)
    return np.where(rr >= r[-1], np.maximum(y[-1], 0.0), np.maximum(out, 0.0))


def add_entropy_core(p: Params, ks: HaloGas, f_core: float, verbose=True, zero_edge=False,
                     mdm=None) -> HaloGas:
    """Re-solve the halo gas for a different entropy profile or a zero-pressure outer edge.

    Entropy is K = P / rho^(5/3).  The K&S halo has K_ks(r) rising outward.
    ``f_core`` adds a constant floor, K(r) = K_ks(r) + f_core * K_ks(r200), which
    flattens the central density into a core and lowers the central pressure.
    Hydrostatic equilibrium dP/dr = -rho G M_dm/r^2 with P = K rho^(5/3) integrates in
    closed form,

        P(r)^(2/5) = P_out^(2/5) + (2/5) int_r^{r_out} K^(-3/5) G M_dm / r'^2 dr'.

    By default P_out is the free parameter, fixed by requiring the gas mass inside
    r200 to stay fbaryon*M200.  f_core = 0 returns the K&S solution exactly.

    ``zero_edge=True`` instead imposes P = 0 at R_max = rtrunc * r200 (as in
    Donnert 2013, where the gas is hydrostatic only out to R_max and the pressure
    vanishes at the boundary), so the gas ends smoothly instead of in a pressure
    cliff and the temperature falls to zero at the edge.  With P_out fixed the
    mass is matched by an overall entropy scale s (K -> sK gives rho -> s^(-3/2) rho).
    Beyond R_max the density is zero.

    ``mdm`` (callable r -> enclosed dark mass [1e10 Msol]) replaces the analytic K&S halo in the
    hydrostatic integral; ``truncated_gas_halo`` uses it with the tapered N-body mass profile.
    """
    if f_core < 0:
        raise ValueError("f_core must be >= 0")
    if f_core == 0 and not zero_edge:
        return ks
    log = print if verbose else (lambda *a, **k: None)
    r = ks.rtab
    x = r / ks.a
    fc = float(_fnorm(ks.kind, ks.conc))
    g = (GGDT * ks.m_ref * _fnorm(ks.kind, x) / fc if mdm is None else GGDT * mdm(r)) / r ** 2
    rho_ks = ks.rho(r)
    p_ks = (GAMMA - 1) * rho_ks * ks.u(r)
    k_ks = p_ks / rho_ks ** (5.0 / 3.0)
    k200 = float(np.interp(np.log(ks.r200), np.log(r), k_ks))
    k = k_ks + f_core * k200
    lr = np.log(r)
    integrand = k ** (-0.6) * g * r                       # dr = r dln r
    cum_in = integrate.cumulative_simpson(integrand, x=lr, initial=0.0)

    from scipy import optimize
    target = ks.mgas
    floor = 1e-300

    def cumulative_mass(rho):
        return integrate.cumulative_simpson(4 * PI * r ** 3 * rho, x=lr, initial=0.0)

    if zero_edge:
        rmax = p.rtrunc * ks.r200
        if not ks.r200 < rmax <= r[-1]:
            raise ValueError("zero_edge needs r200 < rtrunc*r200 <= table end")
        inside = r < rmax
        tail_unit = np.where(inside, np.interp(np.log(rmax), lr, cum_in) - cum_in, 0.0)
        tail_unit = np.maximum(tail_unit, 0.0)

        def profile(ls):
            sc = np.exp(ls)
            pp = (0.4 * tail_unit * sc ** -0.6) ** 2.5     # K -> sK scales the integral by s^-3/5
            rho = np.where(inside, (pp / (sc * k)) ** 0.6, 0.0)
            return pp, rho

        def mass200(ls):
            return float(np.interp(np.log(ks.r200), lr, cumulative_mass(profile(ls)[1])))

        ls = optimize.brentq(lambda v: np.log(mass200(v) / target), -20.0, 20.0, xtol=1e-12)
        pp, rho = profile(ls)
        what = f"zero-pressure edge at {rmax:.4g} kpc ({rmax / ks.r200:.2f} r200), entropy scale {np.exp(ls):.3g}"
    else:
        def profile(p_out):
            pp = (p_out ** 0.4 + 0.4 * (cum_in[-1] - cum_in)) ** 2.5
            return pp, (pp / k) ** 0.6

        def mass200(p_out):
            return float(np.interp(np.log(ks.r200), lr, cumulative_mass(profile(p_out)[1])))

        lo, hi = np.log(p_ks[-1]) - 15, np.log(p_ks[-1]) + 15
        sol = optimize.brentq(lambda lp: np.log(mass200(np.exp(lp)) / target), lo, hi, xtol=1e-12)
        pp, rho = profile(np.exp(sol))
        what = f"free edge pressure, f_core={f_core:g}"
    cum = cumulative_mass(rho)
    u = np.where(rho > 0, pp / ((GAMMA - 1) * np.maximum(rho, floor)), 0.0)
    rho_fn = lambda rr: _loglin_interp(rr, r, rho)
    u_fn = lambda rr: _loglin_interp(rr, r, u)
    temp0 = float(u[0] * U_TO_TEMP * MU_HALO)
    log(f"Halo re-solved ({what}): central density {rho[0] / ks.rho(r[0]):.3g} x K&S; "
        f"P(1 kpc) / K&S = {np.interp(np.log(1.0), lr, pp) / np.interp(np.log(1.0), lr, p_ks):.3g}")
    return HaloGas(ks.kind, ks.a, ks.m200, ks.r200, ks.conc, float("nan"), float("nan"), temp0,
                   float(rho[0]), ks.mgas, r, cum / target, u, rho_fn, u_fn, ks.m_ref,
                   ks.mdm if mdm is None else mdm, rmax if zero_edge else 0.0)


def sample_gas_halo(p: Params, halo: HaloGas, rng, nhalo=None, recentre=False,
                    velocity_noise=True, verbose=True):
    """Poisson-sample the halo.  Returns pos, vel, ug (float32), mpart.

    The particle mass is fixed so that ``nhalo`` particles lie inside r200
    (as the Fortran).  ``velocity_noise`` reproduces the Fortran habit of
    giving each particle a Gaussian velocity with sigma = sqrt(2 u) per
    component (about the sound speed); turn it off for a cold-velocity
    hydrostatic start.
    """
    log = print if verbose else (lambda *a, **k: None)
    nhalo = p.nhalo if nhalo is None else nhalo
    mpart = halo.mgas / nhalo
    router = p.rtrunc * halo.r200
    fm = float(np.interp(router, halo.rtab, halo.mtab))
    nsample = 1 + int(np.floor(fm * halo.mgas / mpart))
    log(f"Sampling {nsample} particles from a Poisson distribution")

    m_enc = rng.random(nsample) * fm
    r = np.interp(m_enc, halo.mtab, halo.rtab)
    az = 2 * PI * rng.random(nsample)
    if isinstance(halo, DiscHaloGas):       # axisymmetric: z = r cos(theta) with the conditional density
        ct = halo.draw_cos_theta(r, rng)
        rc = np.sqrt(1.0 - ct * ct)
        pos = r[:, None] * np.column_stack((rc * np.cos(az), rc * np.sin(az), ct))
        ug = halo.u_rz(np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2])
    else:
        x1 = 2.0 * rng.random(nsample) - 1.0
        rc = np.sqrt(1.0 - x1 * x1)
        pos = r[:, None] * np.column_stack((x1, rc * np.cos(az), rc * np.sin(az)))
        ug = np.interp(r, halo.rtab, halo.utab)
    vel = np.zeros_like(pos)
    vel += p.v_radial * pos / r[:, None]
    if velocity_noise:
        vel += rng.standard_normal(pos.shape) * np.sqrt(2.0 * ug)[:, None]
    if recentre:  # Fortran recentres on the inner half by mass
        order = np.argsort(r)
        pos = pos - 2 * pos[order[: nsample // 2]].sum(0) / nsample
    return pos.astype(np.float32), vel.astype(np.float32), ug.astype(np.float32), mpart
