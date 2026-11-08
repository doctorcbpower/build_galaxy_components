# realtra

![Realtra banner](./realtra_banner.jpg)

Equilibrium initial conditions for galaxy simulations: a **gas disc**, a hot **gas halo** and spherical
**N-body halo/bulge** components, built consistently in one dark-matter potential, stitched together
into one or more **systems** (e.g. for a merger), and optionally converted to an **AREPO** initial condition with a
background grid. 

**Scope.** The value of this package is in building **satellite galaxies** and other gas-rich systems: a small galaxy with its own
gas disc, hot halo and N-body halo, compounded by pressure balance, placed on an orbit around a host, resolved where gas
overlaps, and written with an AREPO background grid. The N-body halo/bulge sampler is a spherical, Eddington-inversion
reference (isotropic, Osipkov–Merritt or constant β); it is there to be understood and to give a working default, not to compete
with dedicated codes. For a host that needs a live halo responding properly to a disc (flattened, self-consistent, action-based or
made-to-measure), build it with GalIC, AGAMA or similar and load that file as a `[[system]]` with `file = "..."` in a merge TOML
(see [Merging systems](#merging-systems)); the satellites and the gas come from here.

The maintained implementation is the Python package [`galaxy_components/`](galaxy_components/). The original
Fortran codes it grew out of are kept in [`legacy/`](legacy/) (with their own README).

```
params.toml ──► disc     gas disc in hydrostatic equilibrium (analytic potential)   ┐
               ──► halo     hot gas halo, hydrostatic in the same dark halo            ├─► system ─► merge ─► IC file
               ──► nbody    spherical isotropic halo/bulge, Eddington inversion        ┘   (compounded)  (orbit, overlap,
                                                                                                      background grid)
```

## Contents
[Install](#install) · [Quick start](#quick-start) · [Parameter file](#parameter-file) · [Commands](#commands) · [GalIC + gas disc](#galic--gas-disc-workflow) · [Running on a cluster](#running-on-a-cluster-arepo) ·
[What each component models](#what-each-component-models) · [Merging systems](#merging-systems) ·
[Use the modules directly](#use-the-modules-directly) · [Tests](#tests) · [Differences from the Fortran](#differences-from-the-fortran) ·
[Limitations](#limitations) · [References](#references)

## Install

Python ≥ 3.11. From a clone of this repository:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[test]"        # numpy>=2, scipy>=1.12, h5py (+ pytest, f90nml); editable, so edits take effect
pytest -q                       # 161 tests (one needs jupyter), ~2 min
galaxy-components --help        # console script (also: python -m galaxy_components ...)
```
The only required dependencies are `numpy`, `scipy` and `h5py`; `f90nml` is optional and only needed to read legacy Fortran
namelist parameter files (`pip install ".[namelist]"`). `pip install .` (non-editable) works too, and `pip wheel .` builds a wheel containing only the package. After installing, the scripts in
[`examples/`](examples/) run from anywhere (without installing, run them from the repository root with `PYTHONPATH=.`).

Units throughout: **kpc, km/s, 10¹⁰ M☉** (so G = 4.3×10⁴ kpc (km/s)² / 10¹⁰ M☉).

## Quick start

```bash
# a whole galaxy (gas disc + gas halo + N-body halo), with an AREPO background grid
galaxy-components galaxy examples/galaxy_params.toml --outdir out --background-grid --halo-core 0.1 --zero-edge

# a two-galaxy parabolic merger from a TOML description
galaxy-components merge examples/merger.toml --outdir out

# add a gas disc to a GalIC snapshot (notebook and scripts in examples/galic_gas_disc/)
galaxy-components check   examples/galic_gas_disc/galic_m1.toml
galaxy-components add-gas examples/galic_gas_disc/galic_m1.toml --snapshot galic.hdf5 --out ic.hdf5

# or script the pieces yourself
python examples/02_gas_halo.py
```
[`examples/galaxy_params.toml`](examples/galaxy_params.toml) is a commented parameter file for a 10¹² M☉ galaxy.

## Parameter file

Parameters are a **TOML** file (standard library, no extra dependency): one table per group, real lists and booleans, comments.
See [`examples/galaxy_params.toml`](examples/galaxy_params.toml) for a commented template.

```toml
[halo]
m200 = 100.0            # virial mass [1e10 Msol]
[gas_disc]
ndisc = 20000
[nbody]
npart_comp = [20000]    # lists may be shorter than 5; unset entries keep their defaults
form_comp = ["nfw"]
```
Files are validated strictly: an unknown table or field, a value of the wrong type, or a field in the wrong table is an error
that names the file and suggests the right table, so a typo can never be silently ignored. Omitted fields take their defaults.

**Legacy namelists.** A Fortran namelist file (`&HALO ... /`, e.g. the old `parameters.txt`) is still read if it does not end in
`.toml` and `f90nml` is installed (`pip install ".[namelist]"`). Convert it once with

```bash
galaxy-components convert parameters.txt -o params.toml      # writes only the non-default values; --all writes everything
```
The legacy Fortran programs in [`legacy/`](legacy/) still need the namelist form ([`examples/galaxy_params.txt`](examples/galaxy_params.txt)).

| table | fields | meaning |
|---|---|---|
| `[halo]` | `m200`, `v200`, `c200`, `qflat`, `rcore` | virial mass [10¹⁰ M☉] (or `v200`), concentration; log-halo flattening/core |
| `[disc]` | `mdisc`, `mstar`, `fdisc`, `fhole`, `fdref`, `sigma0`, `zdisc` | disc mass (>0: fraction of M200; ≤0: absolute), stellar fraction, scale length as a fraction of the halo scale, inner hole, stellar-disc thickness z₀/R_d (GalIC `DiskHeight`; 0 = razor-thin), … |
| `[bulge]` | `mbulge`, `fbulge` | Hernquist bulge mass and scale |
| `[gas_disc]` | `temp`, `ndisc`, `sampling` | gas-disc temperature [K], particle count (0 = no gas disc), placement `"poisson"` (default) or `"sobol"` (quasi-random: about 3x less density shot noise, same distribution) |
| `[gas_halo]` | `fbaryon`, `nhalo`, `rtrunc`, `form`, `v_radial` | baryon fraction, particles inside r200 (0 = no gas halo), sampling edge / r200 (set by `rtrunc_comp` when `istaper`), `"nfw"` or `"hernquist"` |
| `[black_hole]` | `mbh` | central black-hole mass |
| `[nbody]` | `ncomp`, `npart_comp`, `fmass_comp`, `form_comp`, `rscale_comp`, `rtrunc_comp`, `raniso_comp`, `beta_comp`, `istaper` | N-body components (`"nfw"`, `"hernquist"`, `"burkert"`); lists are per component, at most 5; `raniso_comp`/`beta_comp` set the velocity anisotropy ([below](#what-each-component-models)) |
| `[output]` | `disc_file`, `halo_file`, `nbody_file`, `snap_format` | output names; **only HDF5 (`snap_format = 3`) is written** |
| `[model]` | `halo_type`, `galic`, `halo_core`, `halo_zero_edge`, `halo_feels_disc`, `nrbin`, `nzbin`, `log_zmin`, `log_zmax`, `snap_z` | options that replaced the Fortran `-D` compile flags, plus the newer ones |

Parameters that the Fortran namelist carried but that nothing used or implemented (`sigma200`, `lambda_b`, `fdecay`, `glass_file`,
`ispoisson`, `tabfile`, `istabulated`, `exponential_disc`) are not silently accepted: a TOML file that sets one is an error that says
why, and a legacy namelist gives a warning and drops the value (only if it differs from the old default).

`[model]`: `halo_type = "hernquist" | "nfw" | "logarithmic"` (default `"hernquist"`), `galic = true` (GalIC
conventions: the Hernquist halo has total mass `mhalo` and scale `a = r200/c200 · √(2 f(c))`), `halo_core` (entropy
floor, below), `halo_zero_edge` (zero-pressure edge, below), `nrbin`/`nzbin` (gas-disc mesh, default 150 × 500).

## Commands

`galaxy-components <command> PARAMS [--outdir DIR] [--seed N]` (or `python -m galaxy_components ...`)

| command | builds | notable options |
|---|---|---|
| `disc` | gas disc (+ `vrot_components.txt`, `disc_parameters.txt`) | |
| `nbody` | spherical N-body components (+ `distribution_function.txt`) | `--fix-legacy` |
| `halo` | gas halo | `--halo-core F`, `--zero-edge`, `--legacy`, `--cold` |
| `galaxy` / `system` | disc + gas halo + N-body in one file | all `halo` options, plus `--naive-merge`, `--no-nbody`, `--background-grid`, `--boxsize-factor`, `--nngb`, `--name` |
| `merge CONFIG.toml` | several systems stitched into one file | see [Merging systems](#merging-systems) |
| `check` | nothing: rotation curve, gas scale heights, vertical resolution | `--ndisc N` |
| `add-gas` | an existing snapshot plus the gas disc, in one file | `--snapshot`, `--out`, `--ndisc`, `--sampling`, `--strict`, `--background-grid` |

## GalIC + gas disc workflow

The usual use of this package: a collisionless galaxy (dark halo, bulge, stellar disc) already exists, typically from GalIC, and a
gas disc is added in the analytic potential of those components. [`examples/galic_gas_disc/`](examples/galic_gas_disc/) has the
whole workflow twice, with the same parameter file ([`galic_m1.toml`](examples/galic_gas_disc/galic_m1.toml): GalIC Model-M1 plus a gas disc):

* the notebook [`build_gas_disc_ics.ipynb`](examples/galic_gas_disc/build_gas_disc_ics.ipynb) (needs `pip install ".[notebook]"`) goes
  through it with plots: model and rotation curve, vertical resolution against particle number and temperature, the sampled disc
  against the model, the snapshot checks, and the written IC file;
* the command line does the same, step by step or in one go:

```bash
galaxy-components check   examples/galic_gas_disc/galic_m1.toml            # model, rotation curve, h/H by radius
galaxy-components add-gas examples/galic_gas_disc/galic_m1.toml --snapshot galic.hdf5 --out ic.hdf5 --strict
examples/galic_gas_disc/build.sh galic.hdf5 ic.hdf5                         # both (STRICT=1, NDISC=..., BACKGROUND_GRID=1)
```

`add-gas` first compares the snapshot with the analytic model the gas disc is built in (halo, bulge and stellar-disc mass, the
stellar disc's scale length and thickness measured from its particles, its spin direction and centre), because the gas is in
equilibrium only in *that* potential; `--strict` turns a mismatch into an error. It then samples the gas disc, adds it as
PartType0 and writes one GADGET HDF5 file, optionally with an AREPO background grid. `check` reports the vertical resolution
(see [Vertical resolution](#what-each-component-models)) and how much lighter the model halo is than GalIC's, since the gas disc is
taken out of M200 here but is not part of GalIC's halo. The Python functions are in `galaxy_components.workflow`.

## Running on a cluster (Arepo)

[`hpc/`](hpc/) generates a run directory with two SLURM jobs, one that builds the ICs with a background grid and writes the Arepo parameter
file from them, and one that runs Arepo on them (submitted with a dependency). Three ready-made cases: a hot halo plus background grid, the
same plus a cold gas disc, and a GalIC snapshot plus hot halo plus cold gas disc plus background grid.

```bash
cp hpc/cluster.env.example hpc/ozstar.env       # account, modules, Python environment, Arepo binary, resources
python hpc/setup_run.py --galaxy hpc/cases/hot_halo_disc/galaxy.toml --run hpc/cases/hot_halo_disc/run.toml \
    --cluster hpc/ozstar.env --run-dir $SCRATCH/hot_halo_disc --submit             # --scale 0.002 for a quick test
```
See [`hpc/README.md`](hpc/README.md) for the cases, what to check first, and the caveat that the Arepo parameter file is an untested starting point.

## What each component models

**Gas disc** (`disc`). Isothermal gas in hydrostatic equilibrium in the analytic potential of the halo, optional bulge
and stellar exponential disc. At each radius the vertical structure solves
ρ(z) = ρ₀ exp(−(Φ_bg(R,z) − Φ_bg(R,0) + Φ_g(z))/c_s²), Φ_g'' = 4πGρ, with ρ₀ set so the surface density follows the imposed
exponential. The gas's own gravity enters only through that local vertical term and, for the rotation curve, the Freeman
formula; the rotation curve is v_rot² = v_c²(background) + v_c²(gas) + c_s² dlnρ_mid/dlnR. All rings are solved together
with a vectorised RK4 (≈ 2 s for the default mesh). The stellar-disc potential is a Hankel-transform integral.

**Stellar-disc thickness (`zdisc`).** The stellar disc in the background potential is exponential in R and, with `zdisc` = z₀/R_d > 0,
sech²(z/z₀) in height, as in GalIC (`zdisc = DiskHeight`, e.g. 0.2). This matters for the gas disc: a razor-thin stellar disc pulls
on the gas near the plane with a constant 2πGΣ, whereas a thick one gives a force that falls to zero at the plane, so a thin-disc
model makes the gas disc far too thin. For a GalIC Model-M1 stellar disc (z₀ = 0.2 R_d) the half-mass height of the gas disc is
6× the razor-thin value at 0.5 R_d, 5.6× at R_d, 3.9× at 2 R_d, 2.4× at 3 R_d and 1.3× at 5 R_d, and the rotation curve changes by about 2 %. The potential is
the thin-disc Hankel integral plus the exact plane-parallel slab difference 2πGΣ[z₀ ln(2cosh(z/z₀)) − |z|] plus a smooth remainder
tabulated by convolution (agrees with a brute-force convolution to 0.5 % in φ(z) − φ(0); v_c² equals the numerical derivative of φ to
0.2 %). The default `zdisc = 0` is unchanged.

**Vertical resolution.** The gas disc's vertical structure is solved in the continuum limit; a particle code only reproduces it if the
kernel (SPH, 64 neighbours) or cell is small compared with the half-mass height H. When h/H is large the midplane density is
smoothed away, the pressure force falls short of gravity, and the disc collapses and rebounds. In a periodic-patch SPH test
(isothermal column, fixed vertical force) the collapse of H was ~0 at h/H = 1, ~7 % at 2 and ~19 % at 4. `sample_particles` therefore
warns when h/H > 2 (`gas_disc.HRES_MAX`) for more than 5 % of the gas mass, with the radii affected and the `ndisc` that would
resolve 90 % of the gas (h ∝ N^(-1/3)); `gas_disc.resolution` returns h/H ring by ring. `sampling = "sobol"` reduces the density
shot noise about 3× but not this smoothing bias; neither a zero-pressure top edge nor matching u to the SPH density helped in the same test.

**N-body components** (`nbody`). Spherical, multi-component Eddington-type inversion. It sees only the N-body
components themselves: **a disc is neither sampled nor included in its potential**, so the halo is in equilibrium only in
isolation and is perturbed when a disc is added (a true disc + halo equilibrium needs a made-to-measure code such as GalIC).
The gas disc is a perturbation on it, not a self-consistent partner.

*Velocity anisotropy.* By default each component is isotropic. Two exact anisotropic equilibria are available, set per component
in `[nbody]` (at most one of the two for a component):

| setting | model | β(r) | valid range |
|---|---|---|---|
| `raniso_comp = [r_a/r_s, …]` | Osipkov–Merritt, f = f(Q), Q = E − L²/2r_a² | r²/(r² + r_a²): isotropic centre, radial outskirts | β ≥ 0; for an NFW cusp r_a ≳ 0.5 r_s |
| `beta_comp = [β, …]` | constant anisotropy, f = L^(−2β) g(E) | β at all radii, **radial or tangential** | −½ < β ≤ ½ (β = ½ is the radial limit of an r⁻¹ cusp); a cored profile cannot be radially biased |

Both reuse the existing inversion on an effective density (ρ(1 + r²/r_a²) or r^(2β)ρ; for constant β the Abel kernel exponent
becomes ½ + β) and are sampled exactly, in a local frame along r̂ (the isotropic path is unchanged bit for bit). A distribution function
that would go negative is refused with an explanation (the test is whether dropping the negative part would change the density by
more than 0.1 %). Measured from samples, β(r) follows the target to ≈ 0.01 where there are enough particles, and the radial Jeans
equation holds to ≈ 1 %, i.e. these are real equilibria, not just correct second moments.
The isotropic path now also *warns* if its distribution function is negative where it matters (it used to clip it silently).
`raniso_comp` is in units of r_s = r200/c200; `beta_comp` is a Python-only field (the Fortran namelist reader would reject it).

**Cusps.** The distribution function is tabulated on an energy grid. The Fortran's grid is logarithmic in E, which is too coarse
near the centre of a cusp, where Ψ varies by a fraction of a percent over the inner decades: in an NFW halo it gives a velocity
dispersion of 2.0× the Jeans value at r = 0.00015 r_s, 0.73× at 0.0015 r_s and 0.45× at 0.004 r_s. The default grid here has one energy node at
Ψ(r_i) for every radial table node, so it follows the potential (logarithmic in Ψ(r_min) − E in a cusp, in E in the outer halo) and
gives σ/Jeans = 1.00 ± 0.01 from 0.003 kpc outwards at the same number of points; `egrid="log"` (CLI `nbody --log-egrid`)
restores the Fortran grid, and `egrid` also takes an explicit increasing array. Isotropic components are sampled with the same exact
sampler as the anisotropic ones (β = 0); the Fortran's uniform velocity box with a restart at a new radius after 1000 failed tries
moved particles out of a cusp (29 % of the mass inside 0.5 kpc at 10⁶ particles, in the test halo).

**Gas halo** (`halo`). Komatsu & Seljak (2001) polytropic gas, hydrostatic in the dark halo, with Γ matched at r200. With the
default `consistent` mode the gas halo sits in the **same dark halo the disc and N-body parts use**; `--legacy` reproduces the
Fortran (which used a different scale length and mass normalisation). Two options change the halo so it does not overwhelm
the disc:
* `--halo-core F` adds a constant entropy floor K(r) = K_KS(r) + F·K_KS(r200) (K = P/ρ^{5/3}) and re-solves hydrostatic
  equilibrium keeping the gas mass inside r200 (F = 0 is exactly K&S). F = 0.1 raises the fraction of the default disc surviving
  pressure confinement from 19 % to 61 %; beyond ~3 R_d some truncation of a 10⁴ K disc is unavoidable.
* `--zero-edge` (after Donnert 2013) makes the gas hydrostatic only out to R_max = `rtrunc`·r200 with **P(R_max) = 0**, so the
  density and temperature fall smoothly to zero rather than ending in a pressure cliff. The mass inside r200 is matched by an
  overall entropy scale; the sampled halo then holds ≈ 1.01–1.03 × f_b M200 instead of up to 1.7×.
* **Truncated N-body halo (`istaper = true`).** The N-body halo is tapered beyond `rtrunc_comp`·r200 (the largest over the
  components), and the gas halo ends at the same radius with P = 0, whatever `halo_zero_edge` says. The K&S temperature profile is
  the equilibrium of the *untruncated* halo; cut at the taper radius it leaves a pressure cliff and the outskirts expand into the
  tapered potential. So the gas is re-solved in the actual enclosed mass of the N-body model (K&S entropy shape plus `halo_core`,
  mass f_b M200 inside r200): for the default NFW halo tapered at 1.5 r200 it is 1.3× colder at 0.5 r200, 2.4× at r200 and 5.7×
  at 1.3 r200, and falls to zero at the edge. `rtrunc` is then set by `rtrunc_comp`; giving a different value is an error, and
  so is a taper radius inside r200. `--legacy` skips this, as the Fortran did.
* **Halo in the gravity of the disc (`halo_feels_disc = true`).** The disc is a fixed, prescribed mass: exponential in R (mass
  `mdisc`, scale length `rdisc`), treated as razor-thin. The halo density targets are unchanged and the halo adapts to the deeper
  potential, as in made-to-measure codes but without iterating.
  *Gas*: the gas keeps the equation of state P(h) of the spherical solution (h = ∫dP/ρ), and equilibrium in halo + disc is
  h + Φ = const, so ρ(R, z) = ρ(h(R, z)) with h = h_dm(r) − [Φ_disc(R, z) − Φ_disc(pole of the edge)]. The gas is axisymmetric,
  hydrostatic to 1–3 % in both R and z (tested), zero at the same edge, and renormalised to f_b M200 inside r200. It needs a
  zero-pressure edge (`halo_zero_edge` with `rtrunc > 1`, or `istaper`); the compounding rule uses the true P_h(R, z).
  *N-body*: the disc enters the spherical potential through its enclosed mass (the monopole of its potential), so the halo
  density is exactly as specified and the distribution function is inverted in halo + disc (sampled halo obeys the radial Jeans
  equation with the disc mass to 8 %, tested). This under-represents the disc's in-plane pull and does not flatten the halo.
  How much it matters (default NFW halo tapered at 1.5 r200, K&S entropy shape, disc scale length 0.1 r_s):

  | disc mass / M200 | gas density at 5 kpc | default disc mass surviving pressure confinement |
  |---|---|---|
  | 0.02 | × 1.7 | 34 % → 9 % |
  | 0.04 | × 2.4 | 58 % → 18 % |
  | 0.10 | × 3.9 | 78 % → 28 % |

  For the 0.10 disc the N-body halo's velocity dispersion at 1–3 kpc rises from 181 to 421 km/s. A hot halo that ignored the disc
  was much too thin in the centre; a self-consistent one is so dense that a 10⁴ K disc survives only within a few kpc unless
  the halo entropy is raised (`halo_core`).

| halo (default disc) | disc mass kept | disc outer radius |
|---|---|---|
| K&S, free edge | 19 % | 2.6 kpc |
| K&S, zero edge | 17 % | 2.6 kpc |
| core 0.1, free edge | 61 % | 5.7 kpc |
| core 0.1, zero edge | 46 % | 4.5 kpc |

The Fortran also gives each halo particle a Gaussian velocity with σ = √(2u) per component (about the sound speed, ~250 km/s);
this is kept by default, `--cold` removes it.

**Compounding the disc and the halo.** Built independently, the two gas components overlap and the cold disc ignores the pressure of the
hot halo. `galaxy`/`system` therefore apply a pressure-balance rule: the disc exists where its pressure ρ_d c_s² is at least the
halo's P_h(r). The disc is cut at that height (rings whose midplane pressure is below the halo's vanish) and halo particles
inside the disc region are removed. `--naive-merge` just concatenates, as `merge_comps` did. This is approximate: the halo ignores
the disc's gravity and the disc ignores the halo's weight.

**AREPO background grid** (`--background-grid`; port of AREPO's `ADDBACKGROUNDGRID`, `src/add_backgroundgrid`). An oct-tree is
built on the gas; a massless cell is placed at each leaf-node centre; every particle's mass, thermal energy and momentum are spread
onto the cells within its smoothing length with the cubic-spline kernel (conserved exactly; the smoothing length is the distance to
the `--nngb`-th neighbour, since AREPO recomputes it). The box is `--boxsize-factor` (default 5, clamped to AREPO's 2–15) times the
gas extent, all species are shifted to the box centre, `BoxSize` is written to the header, and the file is double precision with a
`Masses` block (zero-mass cells are expected). Printed AREPO values: `BoxSize`, `ReferenceGasPartMass`; AREPO's `MeanVolume` comes
from its internal `GridSize`, which is not reproduced.

## Merging systems

A **system** is one galaxy (gas disc + gas halo + N-body parts, already compounded). `merge` reads a TOML file
([`examples/merger.toml`](examples/merger.toml)) listing any number of systems, each either **built** from a parameter file or
**loaded** from an existing GADGET HDF5 file (including Fortran output), and writes one stitched file.

```toml
output = "merger_ic.hdf5"
seed = 7
[background_grid]            # one grid on the combined gas
enabled = true
[gas_overlap]                # default mode = "pressure"
mode = "pressure"            # or "none"
k = 32
[orbit]                      # positions/velocities of systems 0 and 1 from a Kepler orbit (x-y plane, CM frame)
systems = [0, 1]
r_init = 700.0
r_peri = 30.0
ecc = 1.0
[[system]]
name = "primary"
params = "galaxy_params.toml"         # build ... (a .toml or a legacy namelist file)
[system.options]                      # kwargs of build_system (halo_core, zero_edge, cold_halo, ...)
halo_core = 0.1
[system.overrides]                    # any Params field
ndisc = 40000
[[system]]
name = "secondary"
file = "other.hdf5"                   # ... or load
euler = [0, 60, 0]                    # z-x-z Euler angles [deg]
position = [300, 0, 0]                # (unless an [orbit] sets them)
velocity = [-50, 0, 0]
type_map = { 1 = 2 }                  # keep this system's halo as PartType2
taper_mass = 0.99                     # old merge_comps options: taper_mass,
reduce_stellar_mass = { type = 2, fraction = 0.5 }   #   reduce_stellar_mass,
excise_types = [3]                                   #   excise_types
```

**Gas interfaces between systems.** Where the gas of two systems overlaps, `merge` applies a pairwise pressure-balance rule in the
spirit of Donnert (2013), who resamples the overlapping gas rather than adding it: each gas particle's pressure
P = (γ−1)ρu is estimated from its `k` nearest neighbours, and a particle is removed if any other system's gas would have a higher
pressure at that point, so the higher-pressure gas wins. It needs only particles, so it works for systems loaded from files; the
price is a noisy estimate and a sharp cut. With `mode = "none"` the gas is superposed and a warning is printed if the gas extents
overlap. Choose `r_init` larger than the sum of the gas extents (`rtrunc·r200` each) for the cleanest start.

## Use the modules directly

The package is a set of independent modules, so any piece can be reused on its own. Self-contained, runnable scripts are in
[`examples/`](examples/) (they need the package installed, or `PYTHONPATH=.` from the repository root); they use
`galaxy_components.presets.example_params()` so no parameter file is needed.

| script | shows |
|---|---|
| [`01_gas_disc.py`](examples/01_gas_disc.py) | potential → vertical structure → rotation curve → particles |
| [`02_gas_halo.py`](examples/02_gas_halo.py) | K&S halo vs entropy core vs zero-pressure edge |
| [`03_nbody_halo.py`](examples/03_nbody_halo.py) | spherical N-body halo and a check against the model mass profile |
| [`04_background_grid.py`](examples/04_background_grid.py) | AREPO background grid for **any** gas particles (`--input` an HDF5 file) |
| [`05_stitch_systems.py`](examples/05_stitch_systems.py) | two galaxies, Kepler orbit, interface resolution, stitching |
| [`06_custom_potential.py`](examples/06_custom_potential.py) | compose potentials; reuse the hydrostatic solver for your own Φ(R,z) |
| [`07_anisotropic_halo.py`](examples/07_anisotropic_halo.py) | radially / tangentially biased halos (Osipkov–Merritt, constant β) and measured β(r) |

```python
import numpy as np
from galaxy_components.presets import example_params
from galaxy_components.system import build_system
from galaxy_components import merge

a = build_system(example_params(ndisc=20000), "a", halo_core=0.1)                  # one galaxy in memory
b = build_system(example_params(m200=30.0, ndisc=8000), "b", halo_core=0.1)
x1, v1, x2, v2 = merge.orbit_two_body(a.meta["m200"], b.meta["m200"], r_init=700, r_peri=30)
placed, _ = merge.resolve_gas_overlap([merge.place(a, x1, v1), merge.place(b, x2, v2, euler=(0, 60, 0))])
species, box = merge.finalise(merge.stitch(placed), background_grid=True)           # one background grid
merge.write_merged("merger.hdf5", species, box)
```

| module | contents |
|---|---|
| `potentials` | `Hernquist`, `NFW`, `LogarithmicHalo`, `ExponentialDisc` with `phi(R,z)`, `vc2(R)` |
| `model` | `build_model(params)` – derived radii/masses and the background components |
| `gas_disc` | `build_gas_disc`, `solve_rings` (any Φ), `sample_particles` |
| `gas_halo` | `build_gas_halo`, `add_entropy_core(..., zero_edge=)`, `sample_gas_halo` |
| `nbody` | `build_nbody`, `sample_nbody` |
| `compound` | `truncate_disc`, `excise_halo` (disc ↔ halo pressure balance) |
| `backgroundgrid` | `add_background_grid` (works on any particle arrays) |
| `system` | `build_system`, `read_system` |
| `workflow` | `model_report`, `compare_snapshot`, `add_gas_to_snapshot`, `standin_snapshot` (the GalIC + gas-disc workflow) |
| `merge` | `place`, `orbit_two_body`, `resolve_gas_overlap`, `stitch`, `finalise`, `merge_from_config` |
| `io`, `params`, `presets` | HDF5 writer, TOML/namelist parameter reader and writer, ready-made parameter sets |

## Tests

```bash
pytest -q              # everything, ~1-2 min
pytest -q -m smoke     # end-to-end runs of every command and example, ~25 s
```
The physics tests check against analytic results: the exponential-disc potential against the Freeman closed form, the vertical solver
against the isothermal sech² sheet, hydrostatic balance of every halo variant, conservation in the background grid, the Kepler orbit's
energy and angular momentum, and the N-body distribution function against a fresh run of the Fortran `mk_nbody_comps` (agreement < 1 %).
The anisotropic halos are checked by measuring β(r) from samples, the radial Jeans equation, the r_a → ∞ and β → 0 limits (which must reproduce the isotropic
f with the same normalisation) and the refusal of unrealisable models.
The smoke tests run each command and example on tiny inputs and check the files they write (consistent header, finite values,
unique IDs, particles inside the box). There is no Fortran oracle for the disc: the Fortran `mk_gas_disc` crashes with an
out-of-bounds read in `fintrp` (`legacy/mk_gas_disc.f90`, `xp(0)`).

## Differences from the Fortran

| area | Fortran | here |
|---|---|---|
| π | 3.14159 | exact |
| compile-time choices (`-DNFW_HALO`, …) | makefile flags | runtime (`[model]`) |
| parameter file | Fortran namelist | TOML (namelists still read; `convert` command) |
| gas halo potential | `a = r200/c200`, mass `M200` | same dark halo as the disc (`--legacy` restores) |
| N-body potential interpolation | extrapolates from the wrong table segment | interpolates correctly |
| N-body taper, several components | taper computed with the global NFW scale, so a component with another scale gets a discontinuous density and a negative f | each component's own scale (identical for one component with `rscale_comp = 1`) |
| negative isotropic f | clipped silently | warning when it changes the density by > 1 % |
| isotropic velocity sampling | uniform velocity box, restart at a new radius after 1000 failures (depletes a cusp) | the exact sampler of the anisotropic path at β = 0 |
| N-body velocity scale | `√(G·mass_comp(1)/r_scale)` | `√(G·mtot/r_scale)` (identical for NFW/Burkert) |
| N-body energy grid | logarithmic in E: σ is 0.45× Jeans at r = 0.004 r_s in a cusp | one node per radial node (`egrid="log"` restores) |
| N-body `drptab` edge fill | mixes in `dptab` differences | same by default; `--fix-legacy` corrects |
| disc vertical positions | snapped to the 500-point z table | interpolated (`snap_z` restores) |
| disc rotation curve | last radius reads past the array (`i.lt.nr` bug) | fixed |
| bulge scale in rotation curve | `fbulge·rs`, inconsistent with its potential | same scale as the potential |
| halo recentring | shifts by the inner half's mean | not applied |
| output | Gadget binary or HDF5 | HDF5 only |

## Limitations

* N-body halos are spherical with no rotation. By default the disc is not in their potential; `halo_feels_disc` adds only its
  monopole (enclosed mass), which under-represents its in-plane pull and does not flatten the halo, so it is a first approximation
  for light discs. A proper disc response is deliberately left to dedicated codes (see Scope). The anisotropy is Osipkov–Merritt or constant β in (−½, ½) only. The isotropic sampler
  measures β ≈ −0.14 in the innermost ~0.05 r_s (about 0.1 % of the particles; probably its radius-redraw rule, not investigated); the anisotropic samplers recover
  their target there to ≈ 0.05.
* Systems are superposed, not re-equilibrated, apart from the gas pressure-balance rules.
* Glass-file sampling and Gadget binary output are not ported; `form_comp='general'` is not supported.
* Tabulated gas-halo profiles and the unused Fortran namelist fields are not supported (see [Parameter file](#parameter-file)).
* The AREPO output has not been run in AREPO.
* A 10⁴ K exponential disc inside a hot halo is always pressure-confined at a few scale lengths; see the table above.

## References

Komatsu & Seljak 2001, MNRAS 327, 1353 (gas halo) · Donnert 2013, arXiv:1311.7066 (finite-boundary gas, interface handling) ·
Hernquist 1990 · Springel, Di Matteo & Hernquist 2005 · Yurin & Springel 2014, MNRAS 444, 62 (GalIC) · Freeman 1970 (exponential
disc) · Weinberger, Springel & Pakmor 2020 (AREPO).
