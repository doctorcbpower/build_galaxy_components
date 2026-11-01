"""The halo (gas and N-body) in the gravity of a fixed exponential disc."""
import numpy as np
import pytest
from scipy import integrate

from galaxy_components import compound, gas_disc, gas_halo, nbody
from galaxy_components.constants import GGDT, GAMMA
from galaxy_components.model import build_model
from galaxy_components.potentials import ExponentialDisc
from galaxy_components.presets import example_params
from galaxy_components.system import build_halo

BASE = dict(form="nfw", halo_type="nfw", mdisc=0.06, mstar=0.9, fdisc=0.1, nhalo=15000, ndisc=2000,
            npart_comp=[150000], fmass_comp=[1.0], rscale_comp=[1.0], ncomp=1, form_comp=["nfw"],
            istaper=True, rtrunc_comp=[1.5])


def test_disc_multipole_matches_the_hankel_integral_and_menc_sphere():
    d = ExponentialDisc(4.0, 2.0)
    for R, z in ((59.0, 1.0), (20.0, 50.0), (1e-3, 61.0), (40.0, 40.0)):
        assert d.phi_multipole(R, z) == pytest.approx(d.phi(R, z)[0, 0], rel=1e-5)
    assert d.menc_sphere(1e4) == pytest.approx(4.0)
    # sphere mass = mass inside R < r of the razor-thin disc
    assert d.menc_sphere(3.0) == pytest.approx(4.0 * (1 - 2.5 * np.exp(-1.5)))
    g = d.phi_grid(np.geomspace(0.05, 500, 40), np.geomspace(0.05, 500, 30))
    assert np.isfinite(g).all() and (np.diff(g[:, 0]) > 0).all()    # smooth, rising outwards


@pytest.fixture(scope="module")
def adapted():
    p = example_params(halo_feels_disc=True, **BASE)
    return p, build_halo(p)[0], build_halo(example_params(halo_feels_disc=False, **BASE))[0]


def test_gas_halo_obeys_2d_hydrostatic_equilibrium_in_halo_plus_disc(adapted):
    p, h, _ = adapted
    d = h.disc
    phid = lambda R, z: d.phi_grid(np.array([R]), np.array([z]))[0, 0]
    P = lambda R, z: h.pressure_rz(np.array([R]), np.array([z]))[0]
    for R, z in ((3, 0.5), (10, 3), (10, 30), (40, 5), (30, 60)):
        r = np.hypot(R, z)
        gdm = GGDT * h.mdm(r) / r ** 2
        dR, dz = 1e-3 * R, 1e-3 * max(z, 1.0)
        FR = -(phid(R + dR, z) - phid(R - dR, z)) / (2 * dR) - gdm * R / r
        Fz = -(phid(R, z + dz) - phid(R, z - dz)) / (2 * dz) - gdm * z / r
        rho = h.rho_rz(np.array([R]), np.array([z]))[0]
        assert (P(R + dR, z) - P(R - dR, z)) / (2 * dR) == pytest.approx(rho * FR, rel=0.03)
        assert (P(R, z + dz) - P(R, z - dz)) / (2 * dz) == pytest.approx(rho * Fz, rel=0.03)


def test_gas_mass_edge_and_compression(adapted):
    p, h, sph = adapted
    assert np.interp(h.r200, h.rtab, h.mtab) * h.mgas == pytest.approx(p.fbaryon * h.m200, rel=1e-3)
    assert h.rho_rz(np.array([1.2 * h.redge]), np.array([0.0]))[0] == 0.0          # nothing beyond the edge
    # the disc pulls the gas in: denser in the plane, flattened (denser at fixed r in the plane than on the axis)
    assert h.rho_rz(np.array([5.0]), np.array([0.0]))[0] > 1.5 * sph.rho(5.0)
    assert h.rho_rz(np.array([5.0]), np.array([0.0]))[0] > 1.1 * h.rho_rz(np.array([1e-3]), np.array([5.0]))[0]


def test_sampled_halo_follows_the_axisymmetric_density(adapted):
    p, h, _ = adapted
    pos, vel, ug, mp = gas_halo.sample_gas_halo(p, h, np.random.default_rng(5), nhalo=60000, verbose=False)
    assert len(pos) > 60000
    r = np.linalg.norm(pos.astype(float), axis=1)
    sel = (r > 5) & (r < 10)
    mu2_sample = np.mean((pos[sel, 2].astype(float) / r[sel]) ** 2)
    rr = np.geomspace(5, 10, 40)[:, None]
    mu = np.linspace(0, 1, 400)[None, :]
    rho = h.rho_rz(rr * np.sqrt(1 - mu ** 2), rr * mu)
    w = rr ** 3 * rho
    mu2_model = (integrate.trapezoid(integrate.trapezoid(w * mu ** 2, mu[0], axis=1), rr[:, 0]) /
                 integrate.trapezoid(integrate.trapezoid(w, mu[0], axis=1), rr[:, 0]))
    assert mu2_sample == pytest.approx(mu2_model, rel=0.03)
    assert mu2_model < 0.33                                                          # flattened (a spherical halo gives 1/3)


def test_adapted_halo_requires_a_zero_pressure_edge():
    p = example_params(halo_feels_disc=True, **dict(BASE, istaper=False, rtrunc_comp=[0.0]))
    with pytest.raises(ValueError, match="zero-pressure edge"):
        build_halo(p)
    build_halo(example_params(halo_feels_disc=True, halo_zero_edge=True, rtrunc=1.5,
                              **dict(BASE, istaper=False, rtrunc_comp=[0.0])))


def test_a_denser_halo_confines_the_disc_more(adapted):
    p, h, sph = adapted
    disc = gas_disc.build_gas_disc(p, build_model(p, verbose=False), verbose=False)
    kept = lambda halo: compound.truncate_disc(disc, halo)[0].mrtab[-1] / disc.mrtab[-1]
    assert kept(h) < kept(sph)


def test_nbody_halo_obeys_the_jeans_equation_with_the_disc_included():
    p = example_params(halo_feels_disc=True, **dict(BASE, istaper=False, rtrunc_comp=[0.0], mdisc=0.1))
    nb = nbody.build_nbody(p, verbose=False)
    c = nbody.sample_nbody(p, nb, np.random.default_rng(2), verbose=False)[0]
    iso = nbody.build_nbody(example_params(halo_feels_disc=False, **dict(BASE, istaper=False, rtrunc_comp=[0.0],
                                                                         mdisc=0.1)), verbose=False)
    assert (nb.mtab > iso.mtab).all()                                     # the disc's mass is in the potential
    r, lr = nb.rtab, np.log(nb.rtab)
    mh = nb.menc_tab[0]
    rho = np.gradient(mh, lr) / (4 * np.pi * r ** 3)
    # isotropic Jeans: rho sigma_r^2 (r) = int_r^inf rho G M_tot / r'^2 dr'
    integrand = rho * GGDT * nb.mtab / r ** 2
    cum = integrate.cumulative_simpson(integrand * r, x=lr, initial=0.0)    # dr = r dlnr
    tail = cum[-1] - cum
    sigr2 = tail / rho
    pos = c["pos"].astype(float)
    vel = c["vel"].astype(float)
    rs = np.linalg.norm(pos, axis=1)
    vr = (pos * vel).sum(1) / rs
    rsc = nb.r200 / p.c200
    edges = rsc * np.array([0.1, 0.2, 0.4, 0.8, 1.6, 3.0])
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (rs > lo) & (rs < hi)
        mid = np.sqrt(lo * hi)
        assert vr[s].var() == pytest.approx(np.interp(np.log(mid), lr, sigr2), rel=0.08)
