"""Radially and tangentially biased N-body halos: Osipkov-Merritt and constant beta.

Samples an isotropic halo and three anisotropic ones from exact distribution functions, then
measures beta(r) = 1 - (s_theta^2 + s_phi^2) / (2 s_r^2) from the particles.

    PYTHONPATH=. python examples/07_anisotropic_halo.py
"""
import numpy as np

from galaxy_components import nbody
from galaxy_components.presets import example_params

models = {
    "isotropic": {},
    "Osipkov-Merritt r_a = 2 r_s": dict(raniso_comp=[2.0, 0, 0, 0, 0]),   # beta = r^2/(r^2 + r_a^2), >= 0
    "constant beta = +0.3 (radial)": dict(beta_comp=[0.3, 0, 0, 0, 0]),    # -1/2 < beta < 1/2
    "constant beta = -0.3 (tangential)": dict(beta_comp=[-0.3, 0, 0, 0, 0]),
}

for label, kw in models.items():
    p = example_params(npart_comp=[200000, 0, 0, 0, 0], **kw)
    nb = nbody.build_nbody(p, verbose=False)
    comp = nbody.sample_nbody(p, nb, np.random.default_rng(7), verbose=False)[0]
    rs = nb.r200 / p.c200
    r, beta, n = nbody.measure_beta(comp["pos"], comp["vel"], rs * np.array([0.2, 0.5, 1, 2, 4, 8]))
    print(f"\n{label}   ({len(comp['pos'])} particles)")
    print("   r / r_s :", "  ".join(f"{x / rs:6.2f}" for x in r))
    print("   beta    :", "  ".join(f"{x:6.3f}" for x in beta))
    if nb.kinds[0] == "om":
        print("   target  :", "  ".join(f"{x ** 2 / (x ** 2 + nb.aniso[0] ** 2):6.3f}" for x in r))

# Not every anisotropy is a valid equilibrium: the distribution function must stay >= 0.
print("\nrefused configurations:")
for label, kw in {"Osipkov-Merritt r_a = 0.1 r_s (too radial for a cusp)": dict(raniso_comp=[0.1, 0, 0, 0, 0]),
                  "cored profile with beta = +0.3": dict(form_comp=["burkert", "", "", "", ""], beta_comp=[0.3, 0, 0, 0, 0])}.items():
    try:
        nbody.build_nbody(example_params(**kw), verbose=False)
    except ValueError as e:
        print(f"  {label}:\n     {str(e)[:150]}...")
