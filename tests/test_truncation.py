"""Truncated N-body halo -> gas halo re-solved in the tapered potential, and loud removed parameters."""
import numpy as np
import pytest

from galaxy_components import nbody
from galaxy_components.constants import GAMMA, GGDT
from galaxy_components.gas_halo import U_TO_TEMP, MU_HALO, build_gas_halo
from galaxy_components.params import read_params
from galaxy_components.presets import example_params
from galaxy_components.system import build_halo, build_system


def tapered(**kw):
    base = dict(form="nfw", nhalo=4000, npart_comp=[2000], fmass_comp=[1.0], rscale_comp=[1.0],
                ncomp=1, form_comp=["nfw"], istaper=True, rtrunc_comp=[1.5])
    base.update(kw)
    return example_params(**base)


def test_gas_halo_is_hydrostatic_in_the_tapered_potential_and_ends_at_the_taper_radius():
    p = tapered()
    nb = nbody.build_nbody(p, verbose=False)
    h, pg = build_halo(p, nb=nb)
    assert pg.rtrunc == 1.5
    sel = (h.rtab > 2.0) & (h.rtab < 1.1 * h.r200)
    r = h.rtab[sel]
    rho, u = h.rho(r), h.u(r)
    pres = (GAMMA - 1) * rho * u
    m = np.exp(np.interp(np.log(r), np.log(nb.rtab), np.log(nb.mtab)))
    grav = rho * GGDT * m / r ** 2
    resid = np.abs(np.gradient(pres, np.log(r)) / r + grav) / grav
    assert resid[1:-1].max() < 0.01
    edge = 1.5 * h.r200
    assert (GAMMA - 1) * h.rho(edge) * h.u(edge) < 1e-4 * np.interp(h.r200, r, pres)
    assert h.rho(1.6 * h.r200) == 0
    assert np.interp(h.r200, h.rtab, h.mtab) * h.mgas == pytest.approx(p.fbaryon * h.m200, rel=1e-4)


def test_truncated_gas_is_colder_than_the_free_edge_halo_it_replaces():
    p = tapered()
    h, pg = build_halo(p)
    ks = build_gas_halo(pg, verbose=False)
    for x in (0.5, 1.0, 1.3):
        assert ks.u(x * h.r200) > 1.2 * h.u(x * h.r200)
    # the free-edge halo ends in a pressure cliff; the re-solved one does not
    assert ks.pressure(1.49 * h.r200) > 1e3 * h.pressure(1.49 * h.r200)


def test_untapered_halo_is_unchanged():
    p = tapered(istaper=False, rtrunc_comp=[0.0])
    h, pg = build_halo(p)
    ref = build_gas_halo(p, verbose=False)
    assert pg is p and h.temp0 == ref.temp0 and h.rho(50.0) == ref.rho(50.0)


def test_gas_truncation_radius_errors_are_explicit():
    with pytest.raises(ValueError, match="conflicts with the N-body taper radius"):
        build_halo(tapered(rtrunc=1.2))
    with pytest.raises(ValueError, match="beyond r200"):
        build_halo(tapered(rtrunc_comp=[0.8]))
    with pytest.raises(ValueError, match="rtrunc_comp > 0"):
        nbody.build_nbody(tapered(rtrunc_comp=[0.0]), verbose=False)


def test_system_gas_halo_stays_inside_the_taper_radius():
    p = tapered(ndisc=0, nhalo=3000)
    s = build_system(p, seed=3)
    r = np.linalg.norm(s.species[0]["pos"], axis=1)
    assert r.max() <= 1.5 * s.meta["r200"] * 1.001
    assert (r > s.meta["r200"]).any()


@pytest.mark.parametrize("field,group,value", [
    ("fdecay", "nbody", "3000.0"), ("lambda_b", "gas_halo", "0.1"), ("sigma200", "halo", "1.0"),
    ("ispoisson", "misc", "false"), ("glass_file", "misc", '"g.gdt"'), ("exponential_disc", "model", "false")])
def test_removed_parameters_fail_loudly_in_toml(tmp_path, field, group, value):
    f = tmp_path / "p.toml"
    f.write_text(f"[{group}]\n{field} = {value}\n")
    with pytest.raises(ValueError, match=f"'{field}'.*no longer supported"):
        read_params(f)


def test_removed_parameters_warn_in_a_legacy_namelist(tmp_path):
    pytest.importorskip("f90nml")
    f = tmp_path / "old.txt"
    f.write_text("&HALO\n m200=100.\n/\n&MISC\n glass_file='g.gdt'\n ispoisson=.true.\n/\n")
    with pytest.warns(UserWarning, match="glass_file.*ignored") as w:
        p = read_params(f)
    assert p.m200 == 100.0 and len(w) == 1       # ispoisson at its default is not reported
