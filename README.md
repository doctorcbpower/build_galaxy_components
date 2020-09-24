# build_galaxy_components

This repository contains a set of four codes that have been developed by [A/Prof Chris Power (ICRAR/UWA)](https://www.icrar.org/people/cpower/). The purpose of each is highlighted below:

* [`mk_gas_disc.exe`](mk_gas_disc.f90) - this code is used to compute the hydrostatic equilibrium gas distribution of an isothermal disc embedded within a gravitational potential, usually a composite of a stellar bulge, stellar disc, and dark matter halo. 

* [`mk_gas_halo.exe`](mk_gas_halo.f90) - this code is used to compute the hydrostatic equilibrium gas distribution in a gaseous halo, based on the formalism set out in Komatsu & Seljak 2001.

* [`mk_nbody_comps.exe`](mk_nbody_comps.f90) - this code is used to generate equilibrium n-body components of galaxies with arbitary numbers of components. (...)

* [`merge_comps.exe`](merge_comps.f90) - this code is used to merge gas components with N-body components.

## Usage

### Compilation

To use these codes, they first need to be compiled. Begin by modifying the makefile included in the repository. 
```
> vim makefile
```
There are a number of available options in the makefile that can be activated based on the user's requirements. The first of these catagorise the type of halo potential within which you would like to generate your galaxy componets:
```
#OPT += -DNFW_HALO                  % adding the potential of an analytic NFW halo
#OPT += -DHERNQUIST_HALO            % adding the potential of a Hernquist halo 
#OPT += -DLOGARITHMIC_HALO          % adding the potential of a logarithmic halo 
```
Only one of these halo options can be activiated at a time. 
```
#OPT += -DHERNQUIST_BULGE           % adding the potential of a Hernquist bulge
#OPT += -DSTELLAR_DISC_EXPONENTIAL  % adding the potential of an exponential stellar disk
```
These flags specify the potentials of other components in the N-body model that can optionally be included. Neither, either or both of these options can be activated at a time.
```
#OPT += -DEXPONENTIAL_DISC          % adding the potential of an exponential gas disk
```
This specifies the potential for the gas disc to be placed within the N-body components. This flag needs to be activated when running `mk_gas_disc.exe`.
```
#OPT += -DGalIC                     % turn on GalIC constants
```
This flag needs to be activated in order to use constants and definitions that are consistent with the GalIC N-body initial conditions generation code [(Yurin & Springel, 2014)](https://ui.adsabs.harvard.edu/abs/2014MNRAS.444...62Y/abstract). For example, in the case in which you would like to merge a GalIC N-body simulation and gas components generated with this code. 
```
#OPT += -DHDF5                      % use HDF5
```
This flag needs to be activated to read and generate HDF5 outputs. If turned on, it is necessary to specify the paths of the include and lib folders of your HDF5 installation within the makefile. 

Once you have modified the makefile for your system, the code can be compiled and the executable generated:
```
> make
```
This will compile all executable programs. If you would like to modify the executable for alternative combinations of components, run 
```
> make clean
```
and then recompile the executables with the modified makefile. 


### Inputs
```
> vim parameters.txt 
```
To specify the details of the potential within which the components are built, first modify the input parameters. The parameter.txt file contains a series of namelists; the format is as follows: 
```
&HALO
m200=1.e0               % Virial mass, M200, in units of 1e10 Msol
c200=10                 % NFW concentration
/

&DISC
mdisc=0.035             % Disc mass as a fraction of M200
fdisc=0.2               % Disc scale length as a fraction of halo scale radius
mstar=0.9               % The fraction of the disc mass that is contained in the stellar component
sigma0=0.0              % ?
fhole=0.0               % ?
fdref=-0.0              % ?
/

&BULGE
mbulge=0.0              % Bulge mass as a fraction of M200
fbulge=0.0              % Bulge scale radius as a fraction of the halo scale radius
/

&GAS_DISC
temp=1.d4               % Temperature of the gas disc in K
ndisc=200000            % Number of particles in the gas disc  
/

&GAS_HALO
fbaryon=0.16            % Baryon fraction, OmegaB/Omega0
nhalo=20000             % Number of particles in the gas halo
lambda_b=0.0            % ?
v_radial=0.0            % ?
/

&BLACK_HOLE
mbh=.0e-5               % Black hole mass as a fraction of M200
/

&NBODY      
ncomp=0                 % ?
/

&OUTPUT
disc_file='gas_disc.gdt' % Output gas disc Gadget file
halo_file='gas_halo.gdt' % Output gas halo Gadget file
snap_format=1            % Format of output file (1 - Gadget binary snap1, 2 - Gadget binary snap2, 3 - HDF5)
/

&MISC
glass_file='./glass_32.gdt' % Input glass file
ispoisson=.false            % Sample from a Poisson distribution?
/
```
For the `&MISC` options, only one of these options should be specified. To generate your equilibrium gas component, either:
* Sample your particles from a Poisson distribution, using the ran3() rountine from NR in F77.
or
* By replicating a cube of glass and deforming it to generate the required mass distribution. A glass file is provided in this repository for this use. 

### Running the codes

Each code for generating components accepts the parameters from the input file parameters.txt and can be operated in a similar manner:
```
> ./mk_gas_disc.exe ./parameters.txt    % to build a gas disc within the potential specified 

> ./mk_gas_halo.exe ./parameters.txt    % to build a gas halo within the potential specified

> ./mk_nbody_comps.exe ./parameters.txt % to build an N-body model with components specified 
```
To merge an N-body model with a gas component:
```
> ./merge_comps.exe -in ./gas_component -in ./Nbody_components -out ./merged_output
```
It  is necessary to list the gas component file as the first in the list. If the input/output files are HDF5 format, this needs to be specified by:
```
> ./merge_comps.exe -in_hdf5 ./gas_component.hdf5 -in_hdf5 ./Nbody_components.hdf5 -out_hdf5 ./merged_output.hdf5
```
These flags can be modified for Gadget inputs (-in_snap1 or -in_snap2) or Swift inputs (-in_swift). Similarly, outputs can be written in any of these formats (-out_snap1, -out_snap2 or -out_swift).

## Notes

The mean molecular weight is assumed to 4/(1+3*0.76)~1.21 appropriate for a gas disc; 
for the hot halo, it's assumed to be ionized.

At the moment, I assume the N-body component is set up with GalIC - so potentials and 
parameters are calculated to be consistent with this.

