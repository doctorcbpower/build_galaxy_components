"""Analytic mass components, replacing the #ifdef blocks of the Fortran code.

Each component provides
    phi(R, z)  potential [(km/s)^2]            (arrays broadcast)
    vc2(R)     midplane circular speed squared [(km/s)^2]
so a model is simply a list of components whose contributions are summed.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy import interpolate, special

from .constants import GGDT, PI


class Component:
    def phi(self, R, z):
        raise NotImplementedError

    def vc2(self, R):
        raise NotImplementedError


@dataclass
class Hernquist(Component):
    """Hernquist sphere.  ``mass`` is the asymptotic total mass."""
    mass: float
    a: float

    def phi(self, R, z):
        return -GGDT * self.mass / (self.a + np.hypot(R, z))

    def menc(self, r):
        x = np.asarray(r) / self.a
        return self.mass * x * x / (1 + x) ** 2

    def vc2(self, R):
        return GGDT * self.menc(R) / R


@dataclass
class NFW(Component):
    """NFW halo with ``mass`` the mass inside r200 = c * rs."""
    mass: float
    rs: float
    c: float

    @property
    def _fofc(self):
        return np.log(1 + self.c) - self.c / (1 + self.c)

    def phi(self, R, z):
        x = np.hypot(R, z) / self.rs
        return -GGDT * self.mass / (self.rs * self._fofc) * np.log1p(x) / x

    def menc(self, r):
        x = np.asarray(r) / self.rs
        return self.mass / self._fofc * (np.log1p(x) - x / (1 + x))

    def vc2(self, R):
        return GGDT * self.menc(R) / R


@dataclass
class LogarithmicHalo(Component):
    v0: float
    rcore: float
    q: float = 1.0

    def phi(self, R, z):
        return 0.5 * self.v0 ** 2 * np.log(self.rcore ** 2 + R ** 2 + (z / self.q) ** 2)

    def vc2(self, R):
        return self.v0 ** 2 * R ** 2 / (self.rcore ** 2 + R ** 2)


def _expdisc_nodes(R, Rd, nglq=16):
    """Gauss-Legendre nodes/weights in k for the Hankel integral at radius R.

    Panels resolve (i) the log-spaced scale 1/Rd .. 1/z and (ii) the J0(kR)
    oscillations, so the same rule is valid for all z.
    """
    kmax = 150.0 / Rd
    geom = np.geomspace(1e-4 / Rd, kmax, 40)
    # J0(kR) has period 2pi/R; 16-point panels two periods wide are converged to ~1e-14
    nlin = int(np.ceil(kmax * R / (4 * PI))) + 1
    lin = np.linspace(0.0, kmax, nlin + 1) if nlin > 1 else np.array([0.0, kmax])
    edges = np.unique(np.concatenate(([0.0], geom, lin)))
    xg, wg = np.polynomial.legendre.leggauss(nglq)
    lo, hi = edges[:-1, None], edges[1:, None]
    k = (0.5 * (hi - lo) * xg + 0.5 * (hi + lo)).ravel()
    w = (0.5 * (hi - lo) * wg).ravel()
    return k, w


@dataclass
class ExponentialDisc(Component):
    """Exponential disc, Sigma = M/(2 pi Rd^2) exp(-R/Rd), razor-thin or with a sech^2 thickness.

    Razor-thin:  phi(R,z) = -G M int_0^inf J0(kR) exp(-k|z|) / (1+(k Rd)^2)^{3/2} dk

    ``z0 > 0`` gives rho = Sigma(R) sech^2(z/z0) / (2 z0) (GalIC's ``DiskHeight`` is z0/Rd).  The
    potential is the thin one plus the exact plane-parallel slab difference
    2 pi G Sigma(R) [z0 ln(2 cosh(z/z0)) - |z|], which removes the cusp at z = 0 and carries the
    whole near-plane effect, plus a small smooth remainder (radial-gradient terms) interpolated from a
    table built by convolving the thin potential with the vertical profile.
    """
    mass: float
    rd: float
    z0: float = 0.0
    _delta: object = field(default=None, init=False, repr=False, compare=False)

    @property
    def sigma0(self):
        return self.mass / (2 * PI * self.rd ** 2)

    def _sigma(self, R):
        return self.sigma0 * np.exp(-np.asarray(R, dtype=float) / self.rd)

    def _slab(self, R, z):
        """Thick-minus-thin potential of a plane-parallel slab with the local surface density."""
        x = np.abs(z) / self.z0
        return 2 * PI * GGDT * self._sigma(R) * (self.z0 * (x + np.log1p(np.exp(-2 * x))) - np.abs(z))

    def _build_delta(self, nR=70, nz=48, nodes=48, rmax=20.0):
        """Smooth remainder  int zeta(z') [S(R, z-z') - S(R, z)] dz',  S = thin phi - 2 pi G Sigma |z|."""
        z0 = self.z0
        xg, wg = np.polynomial.legendre.leggauss(nodes)
        zp = 7.0 * z0 * xg
        w = 7.0 * z0 * wg / (2.0 * z0) / np.cosh(zp / z0) ** 2
        Rg = np.geomspace(1e-3 * self.rd, rmax * self.rd, nR)
        zg = np.concatenate(([0.0], np.geomspace(1e-3 * z0, 14.0 * z0, nz)))
        tab = np.empty((nR, zg.size))
        for i, Ri in enumerate(Rg):
            cusp = 2 * PI * GGDT * self._sigma(Ri)
            zz = np.abs(zg[:, None] - zp[None, :])
            vals = ExponentialDisc.phi(self._thin_copy(), Ri, np.concatenate((zz.ravel(), zg)))[0]
            s_conv = vals[:zz.size].reshape(zz.shape) - cusp * zz
            s_at = vals[zz.size:] - cusp * zg
            tab[i] = s_conv @ w - s_at
        self._zscale = 0.05 * z0
        self._zmax = zg[-1]
        self._rmax = Rg[-1]
        self._delta = interpolate.RectBivariateSpline(np.log(Rg), np.arcsinh(zg / self._zscale), tab, kx=3, ky=3)

    def _thin_copy(self):
        return ExponentialDisc(self.mass, self.rd)

    def _delta_at(self, R, z, dR=0):
        if self._delta is None:
            self._build_delta()
        R = np.asarray(R, dtype=float)
        z = np.abs(np.asarray(z, dtype=float))
        zc = np.minimum(z, self._zmax)
        val = self._delta(np.log(np.clip(R, 1e-3 * self.rd, self._rmax)), np.arcsinh(zc / self._zscale),
                          grid=True, dx=dR)
        # beyond the table the remainder falls as (z0/z)^2 in z and vanishes beyond rmax in R
        fall = np.where(z > self._zmax, (self._zmax / np.maximum(z, self._zmax)) ** 2, 1.0)[None, :]
        return val * fall * (R <= self._rmax)[:, None]

    def phi(self, R, z):
        if self.z0 > 0:
            R = np.atleast_1d(np.asarray(R, dtype=float))
            z = np.abs(np.atleast_1d(np.asarray(z, dtype=float)))
            return (ExponentialDisc.phi(self._thin_copy(), R, z) + self._slab(R[:, None], z[None, :]) +
                    self._delta_at(R, z))
        R = np.atleast_1d(np.asarray(R, dtype=float))
        z = np.abs(np.atleast_1d(np.asarray(z, dtype=float)))
        out = np.empty((R.size, z.size))
        for i, Ri in enumerate(R):
            k, w = _expdisc_nodes(Ri, self.rd)
            amp = w * special.j0(k * Ri) / (1 + (k * self.rd) ** 2) ** 1.5
            out[i] = -GGDT * self.mass * (np.exp(-np.outer(z, k)) @ amp)
        return out

    def menc_sphere(self, r):
        """Mass inside a sphere of radius r.  A razor-thin disc inside r is the disc inside R < r; a thick
        one keeps only the part below |z| < sqrt(r^2 - R^2), a fraction tanh(sqrt(r^2 - R^2)/z0)."""
        r = np.asarray(r, dtype=float)
        x = r / self.rd
        if self.z0 <= 0:
            return self.mass * (1.0 - (1.0 + x) * np.exp(-x))
        th, wt = np.polynomial.legendre.leggauss(64)
        th, wt = 0.25 * PI * (th + 1), 0.25 * PI * wt          # R = r sin(theta), theta in [0, pi/2]
        rr = r[..., None]
        integrand = (np.exp(-rr * np.sin(th) / self.rd) * rr * np.sin(th) * rr * np.cos(th) *
                     np.tanh(rr * np.cos(th) / self.z0))
        return 2 * PI * self.sigma0 * (integrand * wt).sum(-1)

    def phi_multipole(self, R, z):
        """Monopole + quadrupole + hexadecapole potential, exact as r -> infinity.

        Q_l = M <R^l> P_l(0), with <R^2> = 6 Rd^2, <R^4> = 120 Rd^4 for an exponential disc.
        Relative error ~ (Rd/r)^6 beyond the last term, ~1e-5 at r = 30 Rd.
        """
        R = np.asarray(R, dtype=float)
        z = np.asarray(z, dtype=float)
        r = np.hypot(R, z)
        mu = z / r
        p2 = 0.5 * (3 * mu ** 2 - 1)
        p4 = (35 * mu ** 4 - 30 * mu ** 2 + 3) / 8
        # Q_l = M <r^l P_l(cos theta)>: r^2 P2 = z^2 - R^2/2,  r^4 P4 = z^4 - 3 z^2 R^2 + 3 R^4 / 8
        z2, z4 = (PI * self.z0) ** 2 / 12, 7 * (PI * self.z0) ** 4 / 240      # <z^2>, <z^4> of sech^2
        r2, r4 = 6.0 * self.rd ** 2, 120.0 * self.rd ** 4
        q2 = self.mass * (z2 - 0.5 * r2)
        q4 = self.mass * (z4 - 3 * z2 * r2 + 0.375 * r4)
        return -GGDT * (self.mass / r + q2 * p2 / r ** 3 + q4 * p4 / r ** 5)

    def phi_grid(self, R, z, rswitch=30.0):
        """phi on the outer product of R and |z|: the Hankel integral within ``rswitch`` scale
        lengths (where it is cheap and exact) and the multipole expansion beyond."""
        R = np.atleast_1d(np.asarray(R, dtype=float))
        z = np.abs(np.atleast_1d(np.asarray(z, dtype=float)))
        rr = np.hypot(R[:, None], z[None, :])
        out = self.phi_multipole(R[:, None], np.maximum(z[None, :], 1e-12))
        near = rr < rswitch * self.rd
        ir = np.where(near.any(1))[0]
        jz = np.where(near.any(0))[0]
        if ir.size:
            exact = self.phi(R[ir], z[jz])
            sub = out[np.ix_(ir, jz)]
            out[np.ix_(ir, jz)] = np.where(near[np.ix_(ir, jz)], exact, sub)
        return out

    def vc2(self, R):
        """Freeman (1970): vc^2 = 4 pi G Sigma0 Rd y^2 [I0K0 - I1K1], y = R/2Rd."""
        R = np.asarray(R, dtype=float)
        y = 0.5 * R / self.rd
        thin = (4 * PI * GGDT * self.sigma0 * self.rd * y * y *
                (special.i0(y) * special.k0(y) - special.i1(y) * special.k1(y)))
        if self.z0 <= 0:
            return thin
        # v_c^2 = R dphi/dR at z = 0: the slab term 2 pi G Sigma z0 ln 2 falls as exp(-R/Rd), plus the remainder
        Rf = np.atleast_1d(R)
        slab = -Rf * 2 * PI * GGDT * self._sigma(Rf) * self.z0 * np.log(2.0) / self.rd
        if self._delta is None:
            self._build_delta()
        dd = self._delta(np.log(np.clip(Rf, 1e-3 * self.rd, self._rmax)), 0.0, grid=True, dx=1)[:, 0] / Rf
        dd = np.where(Rf <= self._rmax, dd * Rf, 0.0)
        return thin + np.reshape(slab + dd, R.shape)
