# Generating ICs and running Arepo on a cluster

`setup_run.py` turns a galaxy model into a **run directory** with two SLURM jobs: one builds the initial conditions
(with an Arepo background grid) and writes the Arepo parameter file from them; the other runs Arepo on those ICs.
The Arepo job is submitted with `--dependency=afterok`, so it starts only if the ICs were made.

```
hpc/
  setup_run.py            creates the run directory and renders the job scripts
  cluster.env.example     per-cluster settings: account, modules, Arepo binary, resources (copy and fill in)
  cases/
    hot_halo/               1. hot gas halo + dark matter particles + background grid
    hot_halo_disc/          2. the same + a cold gas disc (self-consistent: no stellar disc)
    hot_halo_disc_galic/    3. a GalIC snapshot (dark matter, bulge, stellar disc) + hot halo + cold gas disc + background grid
  templates/              make_ics.sbatch.in, arepo.sbatch.in, submit.sh.in, arepo_param.txt.in
```
Each case is a `galaxy.toml` (the model: the same parameter file `galaxy-components` uses) and a `run.toml` (snapshot, duration,
softenings). Cases 1 and 2 build everything from `galaxy.toml`; case 3 builds the gas and adds it to a snapshot you provide.

## Usage

```bash
# once per machine (OzSTAR, NT, ...): copy and fill in  -- nothing in it is filled in for you
cp hpc/cluster.env.example hpc/ozstar.env       # ACCOUNT, MODULES, ENV_SETUP (the Python environment), AREPO_BIN, resources

# 1. a quick test with 0.2 % of the particles, in the foreground on a login node (no sbatch needed)
python hpc/setup_run.py --galaxy hpc/cases/hot_halo/galaxy.toml --run hpc/cases/hot_halo/run.toml \
    --cluster hpc/ozstar.env --run-dir $SCRATCH/hot_halo_test --scale 0.002
bash $SCRATCH/hot_halo_test/make_ics.sbatch

# 2. the real run: write the directory, then submit both jobs
python hpc/setup_run.py --galaxy hpc/cases/hot_halo_disc/galaxy.toml --run hpc/cases/hot_halo_disc/run.toml \
    --cluster hpc/ozstar.env --run-dir $SCRATCH/hot_halo_disc --submit
```
For case 3 set `[run] snapshot = "/path/to/galic.hdf5"` in a copy of its `run.toml` first (the script stops without it), and use your
GalIC disc scale length for `fdisc` in its `galaxy.toml`.

A run directory contains `galaxy.toml`, `merge.toml` (case 3), `make_ics.sbatch`, `arepo.sbatch`, `submit.sh`, `run_settings.json`,
`logs/`, then `ICs/ic.hdf5` (+ `ic.hdf5.arepo.json`, the BoxSize and reference gas mass the Arepo parameters come from),
`arepo_param.txt` and `output/`. If the Arepo job runs out of time, resubmit `arepo.sbatch`: it restarts from `output/restartfiles`.

## What the three cases model

| case | gas | dark matter / stars | notes |
|---|---|---|---|
| `hot_halo` | hot halo, K&S entropy with a core, in the NFW halo, ending at 1.5 r200 with P = 0 | N-body NFW halo, tapered at 1.5 r200 | background grid fills the box (5 x the gas extent) |
| `hot_halo_disc` | + exponential gas disc at 10⁴ K, 3×10⁹ M☉, scale length 3 kpc, quasi-random (Sobol) sampling | N-body NFW halo of mass 0.997 M200 | halo gas and dark matter are in equilibrium with the disc's gravity (`halo_feels_disc`) |
| `hot_halo_disc_galic` | halo + disc as in `hot_halo_disc`, built in the analytic potential of GalIC Model-M1 | GalIC snapshot | stellar disc thickness (`zdisc`) matches GalIC's `DiskHeight` |

## Things to check before a production run

* **The Arepo parameter file is a starting point.** `templates/arepo_param.txt.in` holds the essentials (box, units, softenings, the
  reference gas mass for refinement, output cadence); I could not test it against an Arepo binary. Parameter names and defaults
  differ between versions and depend on your `Config.sh`, so compare it with the example parameter file of your build, or give
  `[arepo] param_template` your own file using the same `@NAME@` placeholders (`IC_FILE`, `OUTPUT_DIR`, `BOXSIZE`, `TIME_MAX`,
  `TIME_BET_SNAP`, `MAX_DT`, `TIMELIMIT_CPU`, `REF_GAS_MASS`, `SOFT_GAS`, `SOFT_DM`, `SOFT_STARS`, `NWRITE`).
* **Gas disc without a stellar disc (case 2).** All the baryons are gas, so Arepo sees the same model the ICs were built in. The disc
  mass is kept at 3×10⁹ M☉ because a cold, gas-only disc is Toomre-unstable at the centre if it is heavier (`galaxy-components disc`
  writes Q(R) to `disc_parameters.txt`). With a stellar disc *as analytic potential only* (`mstar > 0` in case 2) the gas would not be in
  equilibrium in Arepo, which has no such potential: use case 3, where the stars are particles.
* **Vertical resolution of the disc.** At 10⁴ K the disc is a few tens of pc thick, thinner than the gas kernel or cell for any particle
  number you can afford, which makes it collapse and rebound. `galaxy-components check hpc/cases/<case>/galaxy.toml` prints h/H by radius and
  the particle number that would resolve it (see the main README, *Vertical resolution*), and `make_ics.sbatch` warns in its log.
* **Baryon normalisation in case 3.** In the GalIC convention M200 is the Hernquist total mass of the *collisionless* components. The
  hot gas (`fbaryon` x M200) and the gas disc are then added on top, so the system holds more than M200. `check` reports how much lighter the
  model's halo is than GalIC's.
* **Resources.** The defaults in `cluster.env.example` (8 cores, 64 GB, 4 h for the ICs; 2 nodes x 32 tasks, 24 h for Arepo) are guesses.
  Check `sacct -j JOBID --format=MaxRSS,Elapsed` after the first run and adjust. The background-grid step is the memory-hungry one at
  millions of particles.
* The IC job needs the Python environment that has `galaxy-components` installed (`pip install .`), set through `ENV_SETUP`.
