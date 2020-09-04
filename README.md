build_gas_components: containing mk_disc, mk_halo, merge_comps
==============================================================

* mk_disc.exe - code to compute the hydrostatic equilibrium gas 
distribution in an isothermal disc embedded within a gravitational
potential, usually a composite of a stellar bulge, stellar disc, and 
dark matter halo. 

* mk_halo.exe - code to compute the hydrostatic equilibrium gas
distribution in a gaseous halo, based on the formalism set out
in Komatsu & Seljak 2001.

* merge_comps.exe - code to merge gas components with N-body
components.

Usage
=====

% make mk_disc.exe

% vim parameters.txt       ! Edit the parameters for the run

% ./mk_disc.exe ./parameters.txt

% make merge_comps.exe 

% ./merge_comps.exe -gas ./out.gdt -nbody /path/to/galic_nbody_output

% make mk_halo.exe

% vim parameters.txt       ! Edit the parameters for the run

% ./mk_halo.exe ./parameters.txt

% ./merge_comps.exe -in <file1> -in <file2> ... (-in_hdf5|-out_hdf5|-out_snap1|-out_snap2) -out <file>

Inputs
======

The parameters file contains a series of namelists; the format is as follows

&HALO
m200=1.e0               ! Virial mass, M200, in units of 1e10 Msol
c200=6.5                ! NFW concentration
/

&DISC
mdisc=0.035             ! Disc mass fraction, with respect to M200
fdisc=0.208             ! Disc scale length, as a fraction of the halo scale radius, R200/c200
mstar=0.9
/

&BULGE
mbulge=0.0              ! Bulge mass fraction, with respect to M200              
fbulge=0.0              ! Bulge scale length, as a fraction of the halo scale radius, R200/c200
/

&GAS_DISC
temp=1.d4               ! Temperature of the disc, in K
ndisc=200000            ! Number of particles in gas disc
/

&GAS_HALO
fbaryon=0.16            ! Baryon fraction, OmegaB/Omega0
nhalo=20000             ! Number of particles in gas halo
/

&BLACK_HOLE
mbh=0.0                 ! Black hole mass fraction, with respect to M200
/

&OUTPUT
disc_file='gas_disc.gdt'       ! Output gas disc file
halo_file='gas_halo.gdt'       ! Output gas halo file
snap_format=1                  ! Format of output file - (1/2/3) for SnapFormat=1,2, HDF5
/

&MISC
glass_file='./glass_32.gdt'    ! Input glass file
ispoisson=.true.               ! Sample from a Poisson distribution, or deform the glass
/

The mean molecular weight is assumed to 4/(1+3*0.76)~1.21 appropriate for a gas disc; 
for the hot halo, it's assumed to be ionized.

At the moment, I assume the N-body component is set up with GalIC - so potentials and 
parameters are calculated to be consistent with this.

General Comments
================

The user has the option to either lay down gas particles sampled
from a Poisson distribution, generated using the ran3() routine from
NR in F77, or by replicating a cube of glass and deforming it to 
follow the desired mass distribution.

At present the options for the halo, bulge, and disc mass distributions
are

* NFW and Hernquist halos
* Hernquist bulges
* Exponential or Miyamoto-Nagai stelar discs
* Exponential or Power-Law gas discs
