module gadget_header
  integer(kind=4), parameter :: ntype=6
  integer(kind=4), dimension(ntype) :: np,nall
  real(kind=8), dimension(ntype) :: massarr
  real(kind=8) :: expansion=0.,redshift=1.,time=0.
  integer(kind=4) :: flagsfr=0,flagfeedback=0,flagcooling=0,flagdp=0
  integer(kind=4) :: flagstellarage=0, flagmetals=0
  integer(kind=4) :: flagentropy=0, flagdblprc=0
  integer(kind=4), dimension(ntype) :: highword
  integer(kind=4) ::  NumFiles=1
  real(kind=8) :: BoxSize=0.0
  real(kind=8) :: Omega0=1.0,OmegaLambda=0.0,HubbleParam=1.0         
  character :: unused(256-6*4-6*4-6*8-2*8-4*4-4*8-10*4)
  real(kind=8) :: r8dummy
end module gadget_header

