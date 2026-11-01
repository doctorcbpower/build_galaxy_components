"""Exponential stellar disc with a sech^2 vertical profile (GalIC's DiskHeight)."""
import numpy as np
import pytest

from galaxy_components import gas_disc
from galaxy_components.constants import GGDT, PI
from galaxy_components.model import build_model
from galaxy_components.params import read_params
from galaxy_components.potentials import ExponentialDisc
from galaxy_components.presets import example_params


def brute(thick, R, z, n=400):
    """Reference: convolve the thin potential with the vertical profile, cusp subtracted."""
    thin = ExponentialDisc(thick.mass, thick.rd)
    xg, wg = np.polynomial.legendre.leggauss(n)
    L = 12 * thick.z0
    zp = L * xg
    w = L * wg / (2 * thick.z0) / np.cosh(zp / thick.z0) ** 2
    zz = np.abs(z - zp)
    cusp = 2 * PI * GGDT * thin._sigma(R)
    s = thin.phi(R, zz)[0] - cusp * zz
    return s @ w + cusp * thick.z0 * (abs(z) / thick.z0 + np.log1p(np.exp(-2 * abs(z) / thick.z0)))


@pytest.fixture(scope="module")
def disc():
    return ExponentialDisc(4.0, 4.0, 0.8)


@pytest.mark.parametrize("R,z", [(1.0, 0.02), (1.0, 0.3), (4.0, 0.05), (4.0, 1.0), (10.0, 0.2), (10.0, 3.0), (25.0, 0.5)])
def test_potential_matches_the_brute_force_convolution(disc, R, z):
    a = disc.phi(R, z)[0, 0] - disc.phi(R, 0.0)[0, 0]
    b = brute(disc, R, z) - brute(disc, R, 0.0)
    assert a == pytest.approx(b, rel=5e-3)


def test_thin_limit_and_far_field(disc):
    thin = ExponentialDisc(4.0, 4.0)
    tiny = ExponentialDisc(4.0, 4.0, 1e-3)
    assert tiny.phi(8.0, 2.0)[0, 0] == pytest.approx(thin.phi(8.0, 2.0)[0, 0], rel=1e-4)
    assert disc.phi(60.0, 60.0)[0, 0] == pytest.approx(disc.phi_multipole(60.0, 60.0), rel=2e-4)
    assert disc.phi_grid(np.geomspace(0.1, 500, 30), np.geomspace(0.05, 500, 25)).shape == (30, 25)


def test_vertical_force_follows_tanh_in_the_slab_limit():
    d = ExponentialDisc(4.0, 40.0, 0.4)                     # z0 / Rd = 0.01: locally a slab
    R, h = 40.0, 1e-3
    for z in (0.2, 0.4, 0.8):
        fz = (d.phi(R, z + h)[0, 0] - d.phi(R, z - h)[0, 0]) / (2 * h)       # inward force magnitude
        assert fz == pytest.approx(2 * PI * GGDT * d._sigma(R) * np.tanh(z / d.z0), rel=0.04)


def test_circular_speed_is_the_derivative_of_the_potential_and_lower_than_thin(disc):
    thin = ExponentialDisc(4.0, 4.0)
    R = np.array([0.5, 2.0, 4.0, 8.0, 16.0])
    h = 1e-3
    num = [r * (disc.phi(r * (1 + h), 0.0)[0, 0] - disc.phi(r * (1 - h), 0.0)[0, 0]) / (2 * r * h) for r in R]
    assert disc.vc2(R) == pytest.approx(np.array(num), rel=2e-3)
    assert (disc.vc2(R) < thin.vc2(R)).all()


def test_enclosed_mass_of_a_thick_disc():
    d, thin = ExponentialDisc(4.0, 4.0, 0.8), ExponentialDisc(4.0, 4.0)
    r = np.array([0.5, 2.0, 8.0, 200.0])
    m = d.menc_sphere(r)
    assert m[-1] == pytest.approx(4.0, rel=1e-6)
    assert (m[:3] < thin.menc_sphere(r)[:3]).all() and (np.diff(m) > 0).all()


KW = dict(v200=200.0, m200=0.0, c200=10.0, halo_type="hernquist", mdisc=0.04, mstar=0.875, fdisc=0.1, mbulge=0.05,
          fbulge=0.1, ndisc=1000, nrbin=70, nzbin=250, istaper=False, rtrunc_comp=[0.0])


def test_a_thick_stellar_disc_puffs_up_the_gas_disc_but_barely_changes_its_rotation():
    out = {}
    for z in (0.0, 0.2):
        p = example_params(zdisc=z, **KW)
        out[z] = gas_disc.build_gas_disc(p, build_model(p, verbose=False), verbose=False)
    thin, thick = out[0.0], out[0.2]
    rd = thin.model.rdisc
    for fr, lo in ((0.5, 3.0), (1.0, 3.0), (2.0, 2.0)):
        i = np.argmin(abs(thin.rtab - fr * rd))
        h0 = np.interp(0.5, thin.mztab[i], thin.ztab)
        h1 = np.interp(0.5, thick.mztab[i], thick.ztab)
        assert h1 > lo * h0
    i = np.argmin(abs(thin.rtab - 20 * rd))
    assert np.interp(0.5, thick.mztab[i], thick.ztab) < 1.5 * np.interp(0.5, thin.mztab[i], thin.ztab)
    assert abs(thick.vrot[40] / thin.vrot[40] - 1) < 0.05


def test_zdisc_is_a_toml_parameter(tmp_path):
    f = tmp_path / "p.toml"
    f.write_text("[disc]\nmdisc = 0.04\nzdisc = 0.2\n")
    assert read_params(f).zdisc == 0.2
