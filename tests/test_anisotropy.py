"""Velocity anisotropy of the spherical N-body components (Osipkov-Merritt and constant beta)."""
from pathlib import Path

import numpy as np
import pytest

from galaxy_components import nbody
from galaxy_components.constants import GGDT, PI
from galaxy_components.params import read_params

DATA = Path(__file__).parent / "data"
N = 300_000


def make(n=N, **kw):
    p = read_params(DATA / "nbody_ref_params.txt")
    p.npart_comp = [n, 0, 0, 0, 0]
    for k, v in kw.items():
        setattr(p, k, v)
    return p


@pytest.fixture(scope="module")
def models():
    """Build and sample each model once; share between tests."""
    out = {}
    for name, kw in {"iso": {}, "om": dict(raniso_comp=[2.0, 0, 0, 0, 0]),
                     "tangential": dict(beta_comp=[-0.3, 0, 0, 0, 0]),
                     "radial": dict(beta_comp=[0.3, 0, 0, 0, 0])}.items():
        p = make(**kw)
        nb = nbody.build_nbody(p, verbose=False)
        c = nbody.sample_nbody(p, nb, np.random.default_rng(11), verbose=False)[0]
        out[name] = (p, nb, c)
    return out


def bins(nb, p):
    rs = nb.r200 / p.c200
    return rs, rs * np.logspace(np.log10(0.2), np.log10(8), 11)


@pytest.mark.parametrize("name,target", [
    ("om", lambda r, rs: r ** 2 / (r ** 2 + (2.0 * rs) ** 2)),
    ("tangential", lambda r, rs: np.full_like(r, -0.3)),
    ("radial", lambda r, rs: np.full_like(r, 0.3)),
    ("iso", lambda r, rs: np.zeros_like(r)),
])
def test_measured_beta_matches_target(models, name, target):
    p, nb, c = models[name]
    rs, edges = bins(nb, p)
    r, beta, cnt = nbody.measure_beta(c["pos"], c["vel"], edges)
    good = cnt > 5000
    assert good.sum() >= 6
    assert np.max(np.abs(beta[good] - target(r[good], rs))) < 0.04


@pytest.mark.parametrize("name", ["om", "tangential", "radial"])
def test_density_profile_is_unchanged_by_anisotropy(models, name):
    p, nb, c = models[name]
    r = np.linalg.norm(c["pos"].astype(float), axis=1)
    for rr in (0.1 * nb.r200, 0.4 * nb.r200, nb.r200):
        assert abs(c["mpart"] * (r < rr).sum() / np.interp(rr, nb.rtab, nb.menc_tab[0]) - 1) < 0.02


@pytest.mark.parametrize("name", ["iso", "om", "tangential", "radial"])
def test_sampled_halo_obeys_the_jeans_equation(models, name):
    """d ln(nu s_r^2)/d ln r + 2 beta = -G M(<r)/(r s_r^2): the sample is an equilibrium, not just
    a set with the right second moments."""
    p, nb, c = models[name]
    rs = nb.r200 / p.c200
    pos, vel = c["pos"].astype(float), c["vel"].astype(float)
    r = np.linalg.norm(pos, axis=1)
    vr = (vel * pos).sum(1) / r
    vt2 = (vel ** 2).sum(1) - vr ** 2
    edges = rs * np.logspace(np.log10(0.25), np.log10(8), 13)
    mid = np.sqrt(edges[1:] * edges[:-1])
    nu, s2, beta = np.zeros(12), np.zeros(12), np.zeros(12)
    for b in range(12):
        sel = (r >= edges[b]) & (r < edges[b + 1])
        nu[b] = sel.sum() / (4 * PI / 3 * (edges[b + 1] ** 3 - edges[b] ** 3))
        s2[b] = np.mean(vr[sel] ** 2)
        beta[b] = 1 - 0.5 * np.mean(vt2[sel]) / s2[b]
    lhs = np.gradient(np.log(nu * s2), np.log(mid)) + 2 * beta
    rhs = -GGDT * np.interp(mid, nb.rtab, nb.menc_tab[0]) / (mid * s2)
    inner = slice(2, -1)                              # skip the noisy innermost / edge bins
    assert np.max(np.abs(lhs[inner] / rhs[inner] - 1)) < 0.10     # shot noise ~4-6 % (1 sigma) in the innermost bin at this N


def test_isotropic_output_is_unchanged():
    """Regression guard for the isotropic sampler on both energy grids.  Pinned after the sampler switched from the
    Fortran velocity-box rejection to the exact one (193.857 -> 193.802 km/s, 0.03 %) and the default grid from
    logarithmic in E to one node per radial node (193.802 -> 193.381, 0.2 %: the inner cusp is no longer too hot)."""
    p = read_params(DATA / "nbody_ref_params.txt")
    for egrid, expected in ((None, 193.3808277849816), ("log", 193.80221386399668)):
        nb = nbody.build_nbody(p, verbose=False, egrid=egrid)
        c = nbody.sample_nbody(p, nb, np.random.default_rng(19721), verbose=False)
        assert nb.kinds == ("iso",)
        assert nbody.summarise(nb, c)["sigma_tot"] == pytest.approx(expected, rel=1e-12)


def test_limits_reduce_to_the_isotropic_distribution_function():
    """r_a -> infinity and beta -> 0 give the isotropic f with the *same normalisation* (the second
    checks the Gamma-function prefactor of the constant-beta inversion)."""
    iso = nbody.build_nbody(make(), verbose=False)
    om = nbody.build_nbody(make(raniso_comp=[1e6, 0, 0, 0, 0]), verbose=False)
    be = nbody.build_nbody(make(beta_comp=[1e-6, 0, 0, 0, 0]), verbose=False)
    sel = iso.ftab[0] > 0
    assert sel.sum() > 500
    assert np.allclose(om.ftab[0][sel], iso.ftab[0][sel], rtol=1e-5)
    assert np.allclose(be.ftab[0][sel], iso.ftab[0][sel], rtol=1e-3)       # O(beta) away from the limit
    assert nbody.beta_prefactor(0.0) == pytest.approx(1 / (np.sqrt(8) * PI ** 2), rel=1e-12)


def test_generalised_abel_integral_reduces_to_the_eddington_kernel():
    rng = np.random.default_rng(0)
    psi = np.linspace(0, 1, 200)
    g = rng.random(200)
    e = np.linspace(0.05, 1.0, 30)
    assert np.allclose(nbody._abel_integral(e, psi, g, 0.5), nbody._fabel_integral(e, psi, g), rtol=1e-12)
    # kernel exponent: int_0^E (E-s)^(gam-1) ds = E^gam / gam
    one = nbody._abel_integral(e, psi, np.ones(200), 0.3)
    assert np.allclose(one, e ** 0.3 / 0.3, rtol=1e-12)


@pytest.mark.parametrize("kw,fragment", [
    (dict(raniso_comp=[1.0, 0, 0, 0, 0], beta_comp=[0.2, 0, 0, 0, 0]), "not both"),
    (dict(beta_comp=[0.6, 0, 0, 0, 0]), "-1/2 < beta <= 1/2"),
    (dict(beta_comp=[-0.7, 0, 0, 0, 0]), "-1/2 < beta <= 1/2"),
    (dict(beta_comp=[-0.5, 0, 0, 0, 0]), "-1/2 < beta <= 1/2"),
    (dict(raniso_comp=[-1.0, 0, 0, 0, 0]), ">= 0"),
    (dict(raniso_comp=[0.05, 0, 0, 0, 0]), "not realisable"),         # far too radial for an NFW cusp
    (dict(form_comp=["burkert", "", "", "", ""], beta_comp=[0.3, 0, 0, 0, 0]), "not realisable"),   # a core cannot be radial
])
def test_unphysical_anisotropy_is_refused(kw, fragment):
    with pytest.raises(ValueError, match=fragment):
        nbody.build_nbody(make(n=1000, **kw), verbose=False)


def test_cored_halo_accepts_tangential_but_not_radial_bias():
    cored = dict(form_comp=["burkert", "", "", "", ""])
    nbody.build_nbody(make(n=1000, beta_comp=[-0.3, 0, 0, 0, 0], **cored), verbose=False)
    nbody.build_nbody(make(n=1000, raniso_comp=[1.0, 0, 0, 0, 0], **cored), verbose=False)


def test_mixed_components_each_get_their_own_anisotropy():
    """A tapered two-component model: Osipkov-Merritt NFW halo plus a radially biased Hernquist component."""
    p = make(ncomp=2, npart_comp=[40000, 120000, 0, 0, 0], fmass_comp=[1.0, 0.05, 0, 0, 0],
             form_comp=["nfw", "hernquist", "", "", ""], rscale_comp=[1.0, 0.3, 0, 0, 0],
             rtrunc_comp=[2.0, 1.5, 0, 0, 0], raniso_comp=[2.0, 0, 0, 0, 0], beta_comp=[0, 0.2, 0, 0, 0])
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")                      # no 'negative distribution function' warning
        nb = nbody.build_nbody(p, verbose=False)
    assert nb.kinds == ("om", "beta") and nb.aniso[1] == 0.2
    comps = nbody.sample_nbody(p, nb, np.random.default_rng(1), verbose=False)
    assert len(comps) == 2 and all(np.all(np.isfinite(c["vel"])) for c in comps)
    rs = nb.r200 / p.c200
    _, beta, cnt = nbody.measure_beta(comps[1]["pos"], comps[1]["vel"], rs * np.array([0.5, 1, 2, 4]))
    assert np.all(np.abs(beta[cnt > 8000] - 0.2) < 0.06)      # sample scatter is 0.01-0.02 per bin (3 sigma)


def test_clipping_effect_measures_the_density_change():
    e = np.linspace(0.01, 1.0, 400)
    menc = np.linspace(1e-9, 1.0, 400)               # a component whose mass spans the whole table
    psi = e[::-1].copy()                              # Psi decreases with r
    good = e ** 2
    assert nbody.clipping_effect(e, good, 1.5, menc, psi) == 0.0
    dipped = good.copy()
    dipped[100:120] = -0.5 * good[100:120] ** 0.0 * good.max()      # a large negative dip
    assert nbody.clipping_effect(e, dipped, 1.5, menc, psi) > 0.1
    tiny = good.copy()
    tiny[5] = -1e-12 * good.max()                      # numerical dust is harmless
    assert nbody.clipping_effect(e, tiny, 1.5, menc, psi) < 1e-3


def test_negative_isotropic_distribution_function_warns(monkeypatch):
    """The isotropic path used to clip a negative f silently; it must now say so."""
    import warnings
    monkeypatch.setattr(nbody, "clipping_effect", lambda *a, **k: 0.4)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        nbody.build_nbody(make(n=1000), verbose=False)
    assert any("isotropic distribution function is negative" in str(x.message) for x in w)
    # ... and a well-behaved model stays quiet
    monkeypatch.undo()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        nbody.build_nbody(make(n=1000), verbose=False)


# ---------------------------------------------------------------- cusps: beta = 1/2, mass fill, energy grid
def _nfw(n, **kw):
    p = make(n=n, **kw)
    return p, nbody.build_nbody(p, verbose=False)


def test_beta_one_half_is_accepted_for_a_cusp():
    """beta = 1/2 is the radial limit of an r^-1 cusp; it builds, samples and has the right anisotropy."""
    p, nb = _nfw(300_000, beta_comp=[0.5, 0, 0, 0, 0])
    c = nbody.sample_nbody(p, nb, np.random.default_rng(3), verbose=False)[0]
    rs, edges = bins(nb, p)
    _, beta, cnt = nbody.measure_beta(c["pos"].astype(float), c["vel"].astype(float), edges)
    good = cnt > 5000
    assert good.sum() >= 4 and np.all(np.abs(beta[good] - 0.5) < 0.04)


def test_outermost_table_bin_cannot_crash_the_sampler():
    """About 2e-5 of the radii fall in the last table bin, where Psi is below the first energy-grid point
    (f = 0); they are drawn again, so nothing is sampled beyond the grid and the sampler does not raise."""
    for kw in ({}, dict(beta_comp=[0.5, 0, 0, 0, 0]), dict(raniso_comp=[2.0, 0, 0, 0, 0])):
        p, nb = _nfw(400_000, istaper=False, rtrunc_comp=[0.0, 0, 0, 0, 0], **kw)
        c = nbody.sample_nbody(p, nb, np.random.default_rng(5), verbose=False)[0]
        assert c["radius"].max() <= nb.rtab[-2] * (1 + 1e-9)


def test_isotropic_sampler_fills_the_cusp():
    """The Fortran velocity-box rejection restarted particles at a new radius after 1000 failures, which in a
    cusp (f spans ~12 decades) moved them out of the centre: 29 % of the mass inside 0.5 kpc at 1e6 particles."""
    p, nb = _nfw(300_000)
    c = nbody.sample_nbody(p, nb, np.random.default_rng(7), verbose=False)[0]
    r = np.linalg.norm(c["pos"].astype(float), axis=1)
    expected = np.interp(0.5, nb.rtab, nb.menc_tab[0]) / c["mpart"]
    assert 0.75 < (r < 0.5).sum() / expected < 1.3


def _jeans_ratio(nb, rr, rng, n=100_000):
    r, mh = nb.rtab, nb.menc_tab[0]
    rho = np.gradient(mh, r) / (4 * PI * r ** 2)
    integrand = rho * GGDT * mh / r ** 2
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r))))
    with np.errstate(invalid="ignore", divide="ignore"):
        sigr2 = (cum[-1] - cum) / rho                                   # isotropic Jeans: rho sigma_r^2 = int_r^inf rho G M / r'^2
    phi = np.full(n, np.interp(rr, r, nb.ptab))
    vr, vt = nbody._velocities_beta(nb, 0, phi, 0.0, rng)
    return np.sqrt(np.mean(vr ** 2 + vt ** 2) / 3 * nb.vscale ** 2 / np.interp(rr, r, sigr2))


def test_default_energy_grid_resolves_the_cusp():
    """One energy node per radial node follows the potential, so the dispersion of an NFW cusp matches Jeans
    down to 0.003 kpc; the Fortran grid (logarithmic in E) is 2x too hot there and 0.45x too cold at 0.08 kpc."""
    p = make(n=1000)
    nb = nbody.build_nbody(p, verbose=False)
    rng = np.random.default_rng(2)
    assert nb.etab.size == nbody.NTAB - 1 and nb.ftab.shape[1] == nb.etab.size
    for rr in (0.003, 0.01, 0.03, 0.08, 0.3, 2.0):
        assert _jeans_ratio(nb, rr, rng) == pytest.approx(1.0, abs=0.04)
    log = nbody.build_nbody(p, verbose=False, egrid="log")
    assert log.etab.size == nbody.NTAB and _jeans_ratio(log, 0.08, rng) < 0.6      # the known Fortran-grid deficit


def test_default_energy_grid_resolves_anisotropic_cusps():
    rng = np.random.default_rng(4)
    for kw, beta in ((dict(beta_comp=[0.3, 0, 0, 0, 0]), 0.3), (dict(beta_comp=[-0.3, 0, 0, 0, 0]), -0.3)):
        p = make(n=1000, **kw)
        nb = nbody.build_nbody(p, verbose=False)
        phi = np.full(100_000, np.interp(0.03, nb.rtab, nb.ptab))
        vr, vt = nbody._velocities_beta(nb, 0, phi, beta, rng)
        assert 1.0 - 0.5 * np.mean(vt ** 2) / np.mean(vr ** 2) == pytest.approx(beta, abs=0.02)   # beta at 0.03 kpc


def test_energy_grid_argument_is_validated():
    p = make(n=1000)
    with pytest.raises(ValueError, match="egrid"):
        nbody.build_nbody(p, verbose=False, egrid=np.linspace(1.0, 2.0, 50)[::-1])          # not increasing
    with pytest.raises(ValueError, match="egrid"):
        nbody.build_nbody(p, verbose=False, egrid="cusp")
    with pytest.raises(ValueError, match="egrid"):
        nbody.build_nbody(p, verbose=False, egrid=np.linspace(1e9, 2e9, 50))                # outside the potential
