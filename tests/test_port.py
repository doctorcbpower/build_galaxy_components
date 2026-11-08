from pathlib import Path

import numpy as np
from scipy import special

from galaxy_components import gas_disc, nbody
from galaxy_components.constants import GGDT, PI
from galaxy_components.params import read_params
from galaxy_components.potentials import ExponentialDisc

DATA = Path(__file__).parent / "data"


def test_exponential_disc_midplane_potential_matches_freeman():
    d = ExponentialDisc(mass=1.0, rd=2.0)
    R = np.array([0.1, 1.0, 4.0, 20.0, 150.0])
    phi = d.phi(R, np.array([1e-12]))[:, 0]
    y = R / (2 * d.rd)
    exact = -PI * GGDT * d.sigma0 * R * (special.i0(y) * special.k1(y) - special.i1(y) * special.k0(y))
    assert np.allclose(phi, exact, rtol=2e-4)


def test_exponential_disc_vc2_matches_potential_gradient():
    d = ExponentialDisc(mass=1.0, rd=2.0)
    R = np.array([1.0, 3.0, 8.0])
    h = 1e-3
    dphi = (d.phi(R + h, np.array([1e-12])) - d.phi(R - h, np.array([1e-12])))[:, 0] / (2 * h)
    assert np.allclose(R * dphi, d.vc2(R), rtol=1e-3)


def test_isolated_isothermal_sheet_is_sech2():
    # no background: rho = rho0 sech^2(z/z0), z0 = cs^2 / (pi G Sigma) = Sigma/(2 rho0)
    cs, sigma = 9.0, 0.01
    z = 10 ** np.linspace(-6, 1.5, 400)
    rho0, f, cum, _ = gas_disc.solve_ring(sigma, z, np.zeros_like(z), cs)
    z0 = cs ** 2 / (PI * GGDT * sigma)
    assert abs(sigma / (2 * rho0) / z0 - 1) < 2e-3
    assert np.allclose(f, 1 / np.cosh(z / z0) ** 2, rtol=5e-3, atol=1e-6)


def test_nbody_distribution_function_matches_fortran():
    p = read_params(DATA / "nbody_ref_params.txt")
    nb = nbody.build_nbody(p, verbose=False, egrid="log")      # the Fortran oracle uses its own energy grid
    ref = np.loadtxt(DATA / "nbody_ref_distribution_function.txt")
    n = len(ref)
    assert np.allclose(nb.etab[:n], ref[:, 0], rtol=1e-4)
    py = nb.ftab[0][:n] * nb.rscale_comp[0] ** 3
    m = ref[:, 1] > 0
    assert ((py > 0) == m).all()
    assert np.abs(py[m] / ref[m, 1] - 1).max() < 0.015  # ~0.9% observed


def test_vectorised_ring_solver_matches_scipy_reference():
    cs = 9.0
    z = 10 ** np.linspace(-15, 1.3, 500)
    # background: a mix of mild and strong vertical potentials
    R = np.array([0.5, 2.0, 8.0])
    phibz = np.stack([0.5 * k ** 2 * z ** 2 / (1 + z) for k in (30.0, 80.0, 150.0)])
    sig = np.array([0.02, 0.005, 0.0008])
    rho0, f, cum, _ = gas_disc.solve_rings(sig, z, phibz, cs)
    for i in range(3):
        r_ref, f_ref, cum_ref, _ = gas_disc.solve_ring(sig[i], z, phibz[i], cs)
        assert abs(rho0[i] / r_ref - 1) < 2e-4
        assert np.allclose(cum[i, -1], cum_ref[-1], rtol=2e-4)


def test_vectorised_solver_sech2():
    cs, sigma = 9.0, 0.01
    z = 10 ** np.linspace(-6, 1.5, 400)
    rho0, f, cum, _ = gas_disc.solve_rings(np.array([sigma]), z, np.zeros((1, z.size)), cs)
    z0 = cs ** 2 / (PI * GGDT * sigma)
    assert abs(sigma / (2 * rho0[0]) / z0 - 1) < 2e-3


# ---------------------------------------------------------------- gas halo
import pytest

from galaxy_components import backgroundgrid, compound, gas_halo
from galaxy_components.model import build_model


@pytest.mark.parametrize("consistent", [False, True])
def test_gas_halo_is_hydrostatic(consistent):
    p = read_params(DATA / "galaxy_params.txt")
    h = gas_halo.build_gas_halo(p, consistent=consistent, verbose=False)
    r = h.r200 * np.array([0.01, 0.05, 0.2, 0.6, 1.0])
    eps = 1e-4
    dP = (h.pressure(r * (1 + eps)) - h.pressure(r * (1 - eps))) / (2 * eps * r)
    x = r / h.a
    fc = gas_halo._fnorm(h.kind, h.conc)
    # K&S normalisation: M(<r) = m_ref f(x)/f(c), with m_ref = M_halo(<r200)
    m_ref = h.m200 if not consistent else h.m200 * float(fc)
    m_enc = m_ref * gas_halo._fnorm(h.kind, x) / fc
    assert np.allclose(dP / h.rho(r), -GGDT * m_enc / r ** 2, rtol=2e-3)


def test_gas_halo_matches_fortran_reference():
    p = read_params(DATA / "galaxy_params.txt")
    h = gas_halo.build_gas_halo(p, consistent=False, verbose=False)
    assert abs(h.gamma_poly - 1.185) < 1e-9          # Fortran log: 1.185
    assert abs(h.temp0 / 3793166.51 - 1) < 2e-3      # Fortran log: 3793166.51 K
    assert abs(h.mgas - 0.16 * 100.0) < 1e-12


def test_background_grid_conserves_and_tiles():
    rng = np.random.default_rng(3)
    pos = rng.normal(size=(4000, 3)) * np.array([1.0, 1.0, 0.2])
    vel = rng.normal(size=(4000, 3))
    mass = rng.random(4000) + 0.5
    u = rng.random(4000) + 1.0
    out = backgroundgrid.add_background_grid(pos, vel, mass, u, verbose=False)
    m = out["mass"]
    assert np.isclose(m.sum(), mass.sum(), rtol=1e-12)
    assert np.allclose((m[:, None] * out["vel"]).sum(0), (mass[:, None] * vel).sum(0), atol=1e-9)
    assert np.isclose((m * out["u"]).sum(), (mass * u).sum(), rtol=1e-12)
    # leaf nodes are disjoint cubes
    assert out["volume"].sum() < out["boxsize"] ** 3
    lo = out["pos"] - 0.5 * out["volume"][:, None] ** (1 / 3)
    hi = out["pos"] + 0.5 * out["volume"][:, None] ** (1 / 3)
    # no cell centre lies inside another cell
    from scipy.spatial import cKDTree
    t = cKDTree(out["pos"])
    side = out["volume"] ** (1 / 3)
    for i in rng.choice(len(side), 200, replace=False):
        for j in t.query_ball_point(out["pos"][i], r=side[i] * 1.8):
            if j != i:
                assert not np.all((out["pos"][j] > lo[i]) & (out["pos"][j] < hi[i]))


def test_spline_kernel_is_normalised():
    h = 1.7
    r = np.linspace(0, h, 20001)
    assert abs(np.trapezoid(4 * PI * r ** 2 * backgroundgrid.kernel(r / h, h), r) - 1) < 1e-6


def test_disc_truncation_follows_halo_pressure():
    p = read_params(DATA / "galaxy_params.txt")
    p.ndisc = 100
    disc = gas_disc.build_gas_disc(p, build_model(p, verbose=False), verbose=False)
    cold = gas_halo.build_gas_halo(p, verbose=False)
    cold.rho = lambda r: 1e-12 * np.ones_like(np.asarray(r, dtype=float))   # negligible halo
    t, mfrac, _ = compound.truncate_disc(disc, cold)
    assert t.mrtab[-1] > 0.97 * disc.mrtab[-1]
    hot = gas_halo.build_gas_halo(p, verbose=False)
    t2, mfrac2, _ = compound.truncate_disc(disc, hot)
    assert t2.mrtab[-1] < 0.5 * disc.mrtab[-1]            # real halo crushes the outer disc


def test_entropy_core_halo():
    p = read_params(DATA / "galaxy_params.txt")
    p.ndisc = 100
    ks = gas_halo.build_gas_halo(p, verbose=False)
    r = ks.r200 * np.array([0.003, 0.03, 0.3, 1.0])
    near = gas_halo.add_entropy_core(p, ks, 1e-6, verbose=False)
    assert np.allclose(near.rho(r), ks.rho(r), rtol=2e-3)         # f -> 0 recovers K&S
    h = gas_halo.add_entropy_core(p, ks, 0.3, verbose=False)
    assert abs(np.interp(ks.r200, h.rtab, h.mtab) - 1) < 1e-3     # gas mass inside r200 unchanged
    assert h.rho(0.01) < 0.05 * ks.rho(0.01)                       # a real core
    # still hydrostatic in the same dark-matter potential
    rr = ks.r200 * np.array([0.02, 0.1, 0.5])
    eps = 1e-4
    dP = (h.pressure(rr * (1 + eps)) - h.pressure(rr * (1 - eps))) / (2 * eps * rr)
    fc = gas_halo._fnorm(ks.kind, ks.conc)
    m_enc = ks.m_ref * gas_halo._fnorm(ks.kind, rr / ks.a) / fc
    # profile is tabulated (log-linear interpolation), so slopes are good to ~1%
    assert np.allclose(dP / h.rho(rr), -GGDT * m_enc / rr ** 2, rtol=1.5e-2)
    # a gentler halo crushes the disc less
    disc = gas_disc.build_gas_disc(p, build_model(p, verbose=False), verbose=False)
    kept = [compound.truncate_disc(disc, hh)[0].mrtab[-1] for hh in (ks, h)]
    assert kept[1] > 1.5 * kept[0]


def test_zero_pressure_edge_halo():
    p = read_params(DATA / "galaxy_params.txt")
    ks = gas_halo.build_gas_halo(p, verbose=False)
    rmax = p.rtrunc * ks.r200
    for f in (0.0, 0.1):
        h = gas_halo.add_entropy_core(p, ks, f, verbose=False, zero_edge=True)
        assert abs(np.interp(ks.r200, h.rtab, h.mtab) - 1) < 1e-4          # gas mass inside r200 kept
        r = rmax * np.array([0.3, 0.6, 0.9, 0.98, 1.0, 1.05])
        P, T = h.pressure(r), h.u(r) * gas_halo.U_TO_TEMP * gas_halo.MU_HALO
        assert np.all(np.diff(P) <= 0) and P[5] == 0 and P[4] < 1e-6 * P[0]   # monotone, ~zero at the edge (table-cell limited)
        assert np.all(np.diff(T) <= 0) and T[5] == 0 and T[4] < 1e-3 * T[0]  # cold edge
        # hydrostatic away from the edge, in the same dark potential
        rr = ks.r200 * np.array([0.05, 0.3, 0.8])
        eps = 1e-4
        dP = (h.pressure(rr * (1 + eps)) - h.pressure(rr * (1 - eps))) / (2 * eps * rr)
        fc = gas_halo._fnorm(ks.kind, ks.conc)
        m_enc = ks.m_ref * gas_halo._fnorm(ks.kind, rr / ks.a) / fc
        assert np.allclose(dP / h.rho(rr), -GGDT * m_enc / rr ** 2, rtol=2e-2)
        # sampling stops at the edge with no particles beyond it
        pos, vel, ug, mp = gas_halo.sample_gas_halo(p, h, np.random.default_rng(0), nhalo=4000, verbose=False)
        assert np.linalg.norm(pos, axis=1).max() <= rmax * (1 + 1e-6)
        assert np.all(ug >= 0)
