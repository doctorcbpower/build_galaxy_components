module constants
  real(kind=8), parameter :: pi=3.14159
  real(kind=8), parameter :: Mpc_in_m=3.086d22
  real(kind=8), parameter :: Msol_in_kg=1.989d30
  real(kind=8), parameter :: kms_in_ms=1000.
  real(kind=8), parameter :: lunit=1d-3,vunit=1,munit=1.d10
  real(kind=8), parameter :: gmks=6.67d-11
  !  real(kind=8), parameter :: 
  real(kind=8), parameter :: ggdt=gmks*(munit*Msol_in_kg)/(lunit*Mpc_in_m)/(vunit*kms_in_ms)/(vunit*kms_in_ms)
  real(kind=8), parameter :: kb=1.3806d-23  ! Boltzmann constant in MKS
  real(kind=8), parameter :: mp=1.6726d-27  ! Proton mass in kg
  real(kind=8), parameter :: gamma=5./3.
  real(kind=4), parameter :: XH=1.0,YHE=0.,ZM=0.0  ! Abundances of H, He, and metals
  real(kind=8), parameter :: u_to_temp=((gamma-1)*mp/kb)*(kms_in_ms)**2
  real(kind=8), parameter :: joules_to_ergs=1d7
  real(kind=8), parameter :: rhocrit0=27.755d-9 ! Critical density at z=0
  real(kind=8) :: molecular_weight
  real(kind=8), parameter :: deltavir=200.!101.14  
  real(kind=8), parameter :: ln10=dlog(1.d1)
  real(kind=8), parameter :: log10e=dlog10(dexp(1.0d0))

  real(kind=8) :: unit_mass=1000*munit*Msol_in_kg
  real(kind=8) :: unit_length=100*lunit*Mpc_in_m
  real(kind=8) :: unit_time=lunit*Mpc_in_m/vunit/kms_in_ms
  real(kind=8) :: unit_temperature=1.
  real(kind=8) :: unit_current=1.
  integer(kind=4) :: dims=3
end module constants

