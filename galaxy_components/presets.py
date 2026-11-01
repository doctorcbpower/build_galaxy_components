"""Ready-made parameter sets, so the modules can be used without a parameter file."""
from .params import Params


def example_params(**overrides) -> Params:
    """A 1e12 Msol Hernquist/NFW-halo galaxy: gas disc, hot gas halo and an N-body halo.

    Units are kpc, km/s and 1e10 Msol.  Any ``Params`` field can be overridden, e.g.
    ``example_params(ndisc=5000, halo_core=0.1)``.  Coarse disc mesh by default
    (``nrbin``/``nzbin``) so everything runs in seconds.
    """
    p = Params(
        m200=100.0, c200=10.0,                       # halo: M200 [1e10 Msol], concentration
        mdisc=0.04, fdisc=0.076245, mstar=0.9,       # disc: mass/M200, scale length/halo a, stellar fraction
        temp=1e4, ndisc=5000,                        # gas disc: temperature [K], particles
        fbaryon=0.16, nhalo=10000, rtrunc=1.5,       # gas halo: baryon fraction, particles, edge / r200
        form="hernquist",
        ncomp=1, npart_comp=[20000, 0, 0, 0, 0], fmass_comp=[1.0, 0, 0, 0, 0],
        form_comp=["nfw", "", "", "", ""], rscale_comp=[1.0, 0, 0, 0, 0],
        rtrunc_comp=[1.5, 0, 0, 0, 0], istaper=True,
        nrbin=60, nzbin=250,
    )
    for k, v in overrides.items():
        if not hasattr(p, k):
            raise ValueError(f"unknown parameter {k!r}")
        setattr(p, k, v)
    return p
