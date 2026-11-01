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
    assert np.max(np.abs(lhs[inner] / rhs[inner] - 1)) < 0.08     # shot noise ~6 % in the innermost bin at this N


def test_isotropic_output_is_unchanged():
    """Regression guard: the isotropic path must reproduce the pre-anisotropy sample exactly."""
    p = read_params(DATA / "nbody_ref_params.txt")
    nb = nbody.build_nbody(p, verbose=False)
    c = nbody.sample_nbody(p, nb, np.random.default_rng(19721), verbose=False)
    assert nb.kinds == ("iso",)
    assert nbody.summarise(nb, c)["sigma_tot"] == pytest.approx(193.85696362276994, rel=1e-12)


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
    (dict(beta_comp=[0.6, 0, 0, 0, 0]), "-1/2 < beta < 1/2"),
    (dict(beta_comp=[-0.7, 0, 0, 0, 0]), "-1/2 < beta < 1/2"),
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
    assert np.all(np.abs(beta[cnt > 8000] - 0.2) < 0.04)


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
