program main
  use nrutils_modules 
  use gadget_header
  use constants
  use structure
  use test
  use read_params
  use io
  use compute_parameters
  use make_input
  
  implicit none

  integer(kind=4) :: nbox
  
  real(kind=8), allocatable:: lrhobin(:,:)    ! Local density in (r,z) plane
  real(kind=8), allocatable:: rhomid(:)       ! Local density in midplane
  real(kind=8), allocatable :: mrtab(:),mztab(:,:)  ! Tabulated enclosed mass in (r,z)

  real(kind=8), allocatable:: vrot(:)     ! Rotational velocity of gas disc
  real(kind=8), allocatable:: dvrotdr(:)     ! Rotational velocity of gas disc
  real(kind=8) :: sigma_at_r              ! Surface density at radius R
  real(kind=8) :: sigma_ratio             ! Ratio of surface to mid-plane volume densities
  real(kind=8) :: vdisc                   ! Circular velocity due to gas disc
  real(kind=8) :: vcirc                   ! Circular velocity due to spherical components
  real(kind=8) :: vp                      ! Pressure contribution
  logical :: isbh                         ! Is a black hole present?
  real(kind=8) :: mdisc_keep              ! Total mass of disc

  integer(kind=4) :: ndisc      ! Number of particles in disc
  real(kind=8) :: mpart         ! Gas particle mass
  
  integer(kind=4) :: nglass     ! Number of particles in glass file
  integer(kind=4) :: nside     ! Number of particles in glass file

  integer(kind=4) :: nring
  real(kind=8) :: mring,m_enc,rring,rj,r0j,sigma,afac

  real(kind=4), allocatable :: r2(:),xglass(:),yglass(:),zglass(:)
  integer(kind=4), allocatable :: indx(:)

  real(kind=8) :: x,y,z        ! Dummy variables
  real(kind=8) :: lbox         ! Length of box

  real(kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=8) :: xc,yc,zc,zc2,vxc,vyc,vzc,vxc2,vyc2,vzc2
  real(kind=8) :: jx,jy,jz,jtot
  
  integer(kind=4), parameter:: nvert=64

  integer(kind=4) :: nr,nz     ! Dummy index for looping over r,z arrays
  real(kind=8) :: rmin,rmax
  real(kind=8) :: lr,lrmin,lrmax,dlr   ! Log radius
  real(kind=8) :: lz,lzmin,lzmax,dlz   ! Log z
  
  real(kind=8) :: mtot
  real(kind=8) :: dm
  real(kind=8) :: m0,m1,mz,mz0,mz1,z0,z1,mvert
  real(kind=8) :: r0,r1,phi

  real(kind=8) :: rho0old

  integer(kind=4) :: i,j,k,l,m,n,nn
  real(kind=8) :: sum,sum2

  ! Used by odeint
  integer(kind=4) :: ngrid,nok,nbad
  real(kind=8), allocatable ::  xgrid(:),ygrid(:)  
  integer(kind=4), parameter:: nvar=2
  real(kind=8), dimension(nvar) :: vstart
  
  real(kind=8), external :: fintrp

  real(kind=8), external :: dint,fint,phiz

  character(kind=1,len=132) :: paramfile,glassfile,outfile,tabfile
  logical:: fexist
  logical:: ispoisson
  logical:: istabulated  
  integer(kind=4) :: mode           ! mode of output
  external derivs
  
  isbh=.false.

#if defined(HERNQUIST_HALO) && defined(NFW_HALO)
  write(*,*) 'Error: multiple DM halos!'
#endif      

  if(command_argument_count().eq.0) stop 'Usage: mk_disc.exe <parameter_file>'
  
  call get_command_argument(1,paramfile)

  inquire(file=paramfile,exist=fexist)

  if(fexist.eqv..false.) stop 'Error: parameter file does not exist'
  call set_parameters(paramfile,ndisc,outfile,glassfile,istabulated,tabfile,ispoisson,0,mode)
  if(mbh.gt.0.0) isbh=.true.
  
  ! Now compute properties of the disc environment

  if(v200.gt.0.0) then
     mvir=getmvir(v200,deltavir,rhocrit0)
     if(mvir.ne.m200) m200=mvir
  end if
  r200=getrvir(m200,deltavir,rhocrit0)
  rs=r200/c200                                  ! NFW Scale radius

  ! Compute component masses, in units of 1e10 solar masses
  if(mbulge.gt.0) then
     mbulge=mbulge*m200          ! Bulge mass, as fraction of M200
  else
     mbulge=abs(mbulge)          ! Bulge mass, in mass units
  end if

  if(mdisc.gt.0) then
     mdisc=mdisc*m200            ! Total (i.e. gas + stars) disc mass, wrt M200
  else
     mdisc=abs(mdisc)            ! Disc mass, in mass units
  end if
  
  if(mstar.gt.0) then
     mstar=mstar*mdisc           ! Stellar mass, as a fraction of total disc mass
  else
     mstar=abs(mstar)            ! Stellar disc mass, in mass units
  end if
  
#ifndef STELLAR_DISC_EXPONENTIAL
  if(mstar.gt.0) then
     mstar=0.0
     write(*,*) 'Ignoring stellar mass'
  end if
#endif

  if(isbh.eqv..true.) write(*,*) 'Assuming a central black hole...'
  
#ifdef HERNQUIST_HALO
  write(*,*) 'Assuming a Hernquist halo...'
  write(*,*) 'Halo properties:'
  write(*,*) 'M200: ',m200
  write(*,*) 'R200: ',r200
  prof_type=2
#ifndef GalIC
  write(*,*) 'Mtot: ',m200*(1+(1/c200)*sqrt(2*(dlog(1+c200)-c200/(1+c200))))
#else
  write(*,*) 'Mtot: ',m200
#endif
  write(*,*) 'Scale radius:',rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
  write(*,*)
  rdisc=fdisc*rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))       ! in units of kiloparsecs
  rbulge=fbulge*rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))     ! in units of kiloparsecs
  rdisc_hole=fhole*rdisc     ! in units of kpc
#endif
  
#ifdef NFW_HALO
  write(*,*) 'Assuming a NFW halo...'
  write(*,*) 'Halo properties:'
  write(*,*) 'M200: ',m200
  write(*,*) 'R200: ',r200      
  write(*,*) 'Scale radius:',rs
  write(*,*)
  prof_type=1
  rdisc=fdisc*rs                ! in units of kiloparsecs
  rbulge=fbulge*rs              ! in units of kiloparsecs      
  rdisc_hole=fhole*rdisc             ! in units of kpc
#endif            

#ifdef LOGARITHMIC_HALO
  write(*,*) 'Assuming a Logarithmic halo...'
  write(*,*) 'Halo properties:'
  write(*,*) 'M200: ',m200
  write(*,*) 'R200: ',r200  
  write(*,*) 'Core radius:',rcore
  write(*,*) 'q: ',qflat
  write(*,*)
  prof_type=6
  rdisc=abs(fdisc)     
  rbulge=abs(fbulge)
  rdisc_hole=abs(fhole)
  rdisc_ref=abs(fdref)    
#endif
  
#ifdef HERNQUIST_BULGE
  write(*,*) 'Assuming a Hernquist bulge...'
  write(*,*) 'Bulge properties:'
  write(*,*) 'Mtot: ',mbulge
  write(*,*) 'Scale radius:',rbulge
  write(*,*)
#endif  
  
#if defined(STELLAR_DISC_EXPONENTIAL) || defined(STELLAR_DISC_MIYAMOTO_NAGAI)
  write(*,*) 'Assuming an exponential stellar disc...'
  write(*,*) 'Disc properties'
  write(*,*) 'Mstar: ',mstar
  write(*,*) 'Scale radius:',rdisc
#ifdef STELLAR_DISC_MIYAMOTO_NAGAI
  write(*,*) 'Shape parameter:',bdisc
#endif
  write(*,*)
#endif

  mhalo=m200-mbulge-mdisc-mbh
  
  ! Define dimension and bounds of the mesh
  nrbin=150
  lrmin=log10(rdisc)-3.0
  lrmax=log10(rdisc)+log10(30.)
  dlr=(lrmax-lrmin)/real(nrbin-1)
  
  nzbin=500
  lzmin=-15.0
  lzmax=1.3
  dlz=(lzmax-lzmin)/real(nzbin-1)

  allocate(rtab(nrbin))       ! Tabulated radius
  allocate(ztab(nzbin))       ! Tabulated z
  allocate(mrtab(nrbin))      ! Tabulated mass
  allocate(mztab(nrbin,nzbin))  ! Tabulated (normalised) mass in z direction
  
  lr=lrmin
  lz=lzmin
  
  do nr=1,nrbin
     rtab(nr)=10**lr
     lr=lr+dlr
  end do
  
  do nz=1,nzbin
     ztab(nz)=10**lz
     lz=lz+dlz
  end do

  do nr=1,nrbin
     mrtab(nr)=0.0
     do nz=1,nzbin
        mztab(nr,nz)=0.0
     end do
  end do

  write(*,*) 'Radial limits [kpc] :',rtab(1),rtab(nrbin)
  write(*,*) 'Vertical limits [kpc] :',ztab(1),ztab(nzbin)
  write(*,*)
  
  ! Now compute gas disc properties
  write(*,*) 'Rdisc (Gas Disc) [kpc]: ',rdisc
  write(*,*) 'Rhole (Gas Disc) [kpc]: ',rdisc_hole
  write(*,*) 'Rref (Gas Disc) [kpc]: ',rdisc_ref
  rmin=rtab(1)
  if(rdisc_hole.gt.rmin) rmin=rdisc_hole
  rmax=rtab(nrbin)
  if(sigma0.eq.0.0) sigma0=1.0
  call qromb(mdenc,dlog10(rmin),dlog10(rmax),sum)  
  if(sigma0.eq.1.0) sigma0=(mdisc-mstar)/sum
  write(*,*) 'Sigma0 (Gas Disc) [Msol/pc^2]: ',sigma0*1e4,sum*1.4e-3
  write(*,*) 'Gas disc mass [1e10 Msol] :',sigma0*sum  

  call compute_disc_circular_velocity  
  radius=1.0e-3
  rmax=1d3

  call qromb(mdenc,dlog10(rtab(1)),dlog10(rtab(nrbin)),sum)

  molecular_weight=get_mean_molecular_weight(XH,YHE,ZM,temp)
  cs=sqrt(kb*temp/molecular_weight/mp)*1d-3 ! in units of km/s
  write(*,*) 'Sound speed [km/s]: ',cs
  
  ! Now compute total gravitational potential for each of the components
  call compute_background_potential
  
  allocate(lrhobin(nrbin,nzbin),rhomid(nrbin))
  allocate(xgrid(nzbin),ygrid(nzbin))
  
1111 format(a18,i4,4e18.8,a8,i3)

  do nr=1,nrbin
     radius=rtab(nr)
     !     Initial guess at central density         
     call qromb(fint,ztab(1),ztab(nzbin),sum)           
#ifdef EXPONENTIAL_DISC         
#ifndef LOGARITHMIC_HALO
     rho0 = sigma0*exp(-radius/rdisc)/(2.*sum)
#else
     rho0 = sigma0*(radius*radius/(radius*radius+1.))*exp(-(radius-rdisc_ref)/rdisc)/(2.*sum)
#endif
#else
     rho0 = sigma0*(rdisc/radius)/(2.*sum)
#endif         

#ifdef EXPONENTIAL_DISC         
#ifndef LOGARITHMIC_HALO
     sigma_ratio=sigma0*exp(-radius/rdisc)/(2*rho0)
#else
     sigma_ratio=sigma0*(radius*radius/(radius*radius+1.))*exp(-(radius-rdisc_ref)/rdisc)/(2.*rho0)
#endif
#else
     sigma_ratio=sigma0*(rdisc/radius)/(2.*rho0)     
#endif
     write(*,*)     
     write(*,1111) 'Initial guess: ',nr,radius,radius/rdisc,rho0,sigma_ratio
     
     rho0old=0.0
     k=0

     zmax=ztab(nzbin)

     do
        k=k+1
        
        vstart(1)=0.0  ! Boundary conditions in midplane of disc
        vstart(2)=0.0

        call odeint(vstart,nvar,ztab(1),zmax,1.0d-6,&
             &           1.d-16,0.0d0,nok,nbad,derivs,rkqs,xgrid,ygrid,ngrid,&
             &           phiz,fintrp)

        call mqromb(fintrp,ztab(1),zmax,xgrid,ygrid,ngrid,sum,0)
#ifdef EXPONENTIAL_DISC
#ifndef LOGARITHMIC_HALO
        rho0 = sigma0*exp(-radius/rdisc)/(2.*sum)
#else
        rho0 = sigma0*(radius*radius/(radius*radius+1.))*exp(-(radius-rdisc_ref)/rdisc)/(2.*sum)
#endif
#else
        rho0 = sigma0*(rdisc/radius)/(2.*sum)
#endif        
        if(abs(rho0old-rho0).lt.(1.0d-3)) exit

        if(k.gt.300) exit
        
        rho0old = rho0
     end do

#ifdef EXPONENTIAL_DISC         
#ifndef LOGARITHMIC_HALO
     sigma_ratio=sigma0*exp(-radius/rdisc)/(2*rho0)
#else
     sigma_ratio=sigma0*(radius*radius/(radius*radius+1.))*exp(-(radius-rdisc_ref)/rdisc)/(2.*rho0)
#endif
#else
     sigma_ratio=sigma0*(rdisc/radius)/(2.*rho0)     
#endif
     
     write(*,1111) 'Converged result: ',nr,radius,radius/rdisc,rho0,sigma_ratio,', niter: ',k
     write(34,*) nr,radius,radius/rdisc,rho0,sigma_ratio,k
     rhomid(nr)=rho0
     
     do j=1,nzbin
        i=1
        do
           if(xgrid(i).gt.ztab(j).or.i.eq.ngrid) exit
           i=i+1
        end do
        y = (ztab(j)-xgrid(i-1))/(xgrid(i)-xgrid(i-1))*ygrid(i)&
             &           + (xgrid(i)-ztab(j))/(xgrid(i)-xgrid(i-1))*ygrid(i-1)
        lrhobin(nr,j)=log10(rho0*exp(-(phiz(ztab(j))+y)/cs**2))
        call mqromb(fintrp,10**lzmin,ztab(j),xgrid,ygrid,ngrid,sum2,0)
        mztab(nr,j)=sum2/sum
        if(1.0-mztab(nr,j).lt.5.e-3) mztab(nr,j)=1.0
     end do
     mrtab(nr)=(mdisc-mstar)*(1.d0-exp(-rtab(nr)/rdisc)*(1+rtab(nr)/rdisc))
     print *, nr,rtab(nr),mrtab(nr)
  end do

  write(*,*) 'Finished computing density structure...'
  write(*,*)

  ! Compute the circular velocity at the position of each radial bin.  

  allocate(vrot(nrbin))

  open(66,file='vrot_components.txt',status='unknown')
#ifdef LOGARITHMIC_HALO
  write(66,*) '# Logarithmic halo'
#endif  
#ifdef HERNQUIST_HALO
  write(66,*) '# Hernquist halo'
#endif
#ifdef NFW_HALO
  write(66,*) '# NFW halo'
#endif
#ifdef HERNQUIST_BULGE
  write(66,*) '# Hernquist bulge'
#endif
#ifdef STELLAR_DISC_EXPONENTIAL
  write(66,*) '# Exponential disc'
#endif    
  write(66,'(a13,f12.8)') '# rdisc [kpc]:',rdisc
  write(66,'(a19,f12.8)') '# mdisc [1e10 Msol]:',mdisc
  write(66,'(a19,f12.8)') '# mstar [1e10 Msol]:',mstar  
  write(66,'(a27,f12.8)') '# sigma0 [1e10 Msol/kpc^3]:',sigma0
  write(66,'(4a18)') '# Radius [kpc]','Vc [km/s]','Vdisc [km/s]','Vp [km/s]'
  
  do i=1,nrbin
     radius=rtab(i)
     if(i.gt.1.and.i.lt.nr) then
        vp=0.5*cs*cs*(lrhobin(i+1,1)-lrhobin(i-1,1))/dlr
     else if(i.eq.1) then
        vp=cs*cs*(lrhobin(i+1,1)-lrhobin(i,1))/dlr
     else
        vp=cs*cs*(lrhobin(i,1)-lrhobin(i-1,1))/dlr
     end if
     vcirc=0.0

#if defined(HERNQUIST_HALO) || defined(NFW_HALO)
     vcirc=vcirc+ggdt*mhenc(radius)/radius
#endif

#ifdef LOGARITHMIC_HALO
     vcirc=vcirc+circular_velocity_logarithmic_halo(v200,rcore,radius)
#endif     
     
#ifdef HERNQUIST_BULGE         
     vcirc=vcirc+ggdt*mbenc(radius)/radius
#endif
     
     y=0.5*radius/rdisc
     
#ifdef STELLAR_DISC_EXPONENTIAL         
     vcirc=vcirc + 4.*PI*ggdt*(mstar/(mdisc-mstar))*sigma0*rdisc*y*y*(bessi0(y)*bessk0(y)-&
          &        bessi1(y)*bessk1(y))
#endif
     vdisc=4.*PI*ggdt*sigma0*rdisc*y*y*(bessi0(y)*bessk0(y)-bessi1(y)*bessk1(y))
     
     vrot(i)=sqrt(vcirc+vdisc+vp)
     write(66,*) radius, sqrt(vcirc),sqrt(vdisc),sqrt(abs(vp)),vrot(i)
  end do
  close(66)
  write(*,*) 'Finished computing rotational velocity profile...'

  ! Tabulate useful parameters
  
  allocate(dvrotdr(nrbin))

  open(39,file='disc_parameters.txt',status='unknown')
  write(39,'(a13,f12.8)') '# rdisc [kpc]:',rdisc
  write(39,'(a19,f12.8)') '# mdisc [1e10 Msol]:',mdisc
  write(39,'(a19,f12.8)') '# mstar [1e10 Msol]:',mstar  
  write(39,'(a27,f12.8)') '# sigma0 [1e10 Msol/kpc^3]:',sigma0
  write(39,'(4a18)') '# Radius','Sigma(R)','Vrot','Q'
  write(39,'(4a18)') '# [kpc]','[1e10 Msol/kpc^3]','Vrot',''  
  call get_derivative(rtab,vrot,dvrotdr,nrbin)

  do i=1,nrbin
#ifdef EXPONENTIAL_DISC         
     sigma_at_r = sigma0*exp(-rtab(i)/rdisc)
#else
     sigma_at_r = sigma0*(rtab(i)/radius)
#endif         
     write(39,'(5e18.8)') rtab(i),sigma_at_r,vrot(i),&
          & get_qparameter(rtab(i),vrot(i),dvrotdr(i),cs,&
          & sigma0*dexp(-rtab(i)/rdisc)),&
          & sigma_at_r/2/rhomid(i)
  end do

  close(39)

  deallocate(dvrotdr,rhomid)
  
  ! Now impose density profile on glass
  
  ! How extended is the disc? This will set the mass.
  
  lrmax=min(lrmax,log10(rdisc)+1.4)
  
  i=1
  do
     if(rtab(i).ge.10**lrmax) exit
     i=i+1
  end do
  
  mdisc_keep=(rtab(i)-10**lrmax)/(rtab(i)-rtab(i-1))*mrtab(i-1)&
       &     +(10**lrmax-rtab(i-1))/(rtab(i)-rtab(i-1))*mrtab(i)
  
  ! Is there an inner edge to the disc?
  
  r0=fhole*rdisc  ! Inner edge of ring
  
  i=1
  
  do
     if(rtab(i).gt.r0) exit
     i=i+1
  end do
  
  if(i.eq.0) then
     m0=0
  else
     m0=(rtab(i)-r0)/(rtab(i)-rtab(i-1))*mrtab(i-1)&
          &     +(r0-rtab(i-1))/(rtab(i)-rtab(i-1))*mrtab(i)
  end if
  
  ! Given disc mass and number of particles, now know the mass of each particle.
  
  mpart=(mdisc_keep-m0)/real(ndisc)

  ! What is the maximum extent the disc can have, given the mass resolution assumed?
  dm=mdisc_keep-mpart
  
  i=1
  
  do
     if(mrtab(i).gt.dm) exit
     i=i+1
  end do
  
  r1=(mrtab(i)-dm)/(mrtab(i)-mrtab(i-1))*rtab(i-1)&
       &        +(dm-mrtab(i-1))/(mrtab(i)-mrtab(i-1))*rtab(i)

  write(*,*) 'For assumed particle mass, disc radial extent [kpc]:',r1

  if(ispoisson.eqv..false.) then
     inquire(file=glassfile,exist=fexist)
     if(fexist.eqv..false.) stop 'Could not find input glass file...'
     call get_nglass(glassfile,nglass)

     ! First check how many particles are on a side...
     nside=floor(real(nglass)**(1./3.))
     ! ... and then note that you want to keep (nvert/nside) x nglass
     ! in the glass, which you then divide into ndisc to get total
     ! number of particles in slab. Take the square root to get NxN
     ! cubes of glass.
     nbox=1+floor(sqrt(real(ndisc)/(real(nvert)/real(nside)*real(nglass))))

     allocate(xglass(nglass*(nbox+1)*(nbox+1)))
     allocate(yglass(nglass*(nbox+1)*(nbox+1)))
     allocate(zglass(nglass*(nbox+1)*(nbox+1)))
     if(allocated(xglass).eqv..false. .or.&
          & allocated(yglass).eqv..false. .or.&
          & allocated(zglass).eqv..false.)&
          &  stop 'Could not allocate memory!'     
     call get_glass_particles(glassfile,nvert,nbox,nglass,&
          & nglass*(1+nbox)*(1+nbox),xglass,yglass,zglass,1)
  end if
  
  allocate(pos(3,ndisc+1),vel(3,ndisc+1),id(ndisc+1),ug(ndisc))
  
  nn=0

  xc=0.0
  yc=0.0
  zc=0.0
  
  vxc=0.0
  vyc=0.0
  vzc=0.0
  
  vxc2=0.0
  vyc2=0.0
  vzc2=0.0
  
  jx=0.0
  jy=0.0
  jz=0.0
  
  write(*,*) 'Looping over particles...'

  if(ispoisson) then
     iseed=-19723616
     do
        mtot=ran3(iseed)*mdisc_keep   ! Know total mass - random number sets the enclosed mass

        if(nn.gt.ndisc) then
           nn=nn-1
           exit
        end if
        
        i=1
        
        do
           if(mrtab(i).gt.mtot) exit    ! Use this to deduce radius
           i=i+1
        end do
        
        r1=(mrtab(i)-mtot)/(mrtab(i)-mrtab(i-1))*rtab(i-1)&
             &        +(mtot-mrtab(i-1))/(mrtab(i)-mrtab(i-1))*rtab(i)
        
        if(r1.lt.r0) cycle

        nn=nn+1                         ! Increment particle count by 1
        
        vcirc=(mrtab(i)-mtot)/(mrtab(i)-mrtab(i-1))*vrot(i-1)&
             &   +(mtot-mrtab(i-1))/(mrtab(i)-mrtab(i-1))*vrot(i)    ! Compute circular velocity
        
        mvert=ran3(iseed)               ! Now work out corresponding height
        j=1
        
        do
           if(j.ge.nzbin) exit
           if(mztab(i-1,j).gt.mvert) exit
           j=j+1
        end do
        
        z=ztab(j)
        if(1.-2*ran3(iseed).lt.0) z=-1*z            
        phi=2*pi*ran3(iseed)

        pos(1,nn)=r1*cos(phi) 
        pos(2,nn)=r1*sin(phi) 
        pos(3,nn)=z
     
        vel(1,nn)=-vcirc*sin(atan2(pos(2,nn),pos(1,nn)))
        vel(2,nn)=vcirc*cos(atan2(pos(2,nn),pos(1,nn)))
        vel(3,nn)=0.0
        
        id(nn)=nn
        
        ug(nn)=cs*cs/(gamma-1)
     end do
  else
     mtot=m0
     dm=mpart
     do
        mtot=mtot+dm
     
        if(mtot.ge.mdisc_keep) exit

        ! Get the outer radius r1 -- inner radius given by r0
        i=1
        
        do
           if(mrtab(i).ge.mtot) exit
           i=i+1
        end do

        r1=(mrtab(i)-mtot)/(mrtab(i)-mrtab(i-1))*rtab(i-1)&
             &        +(mtot-mrtab(i-1))/(mrtab(i)-mrtab(i-1))*rtab(i)
        
        if(r1.lt.r0) cycle

        nn=nn+1
        
        vcirc=(mrtab(i)-mtot)/(mrtab(i)-mrtab(i-1))*vrot(i-1)&
             &   +(mtot-mrtab(i-1))/(mrtab(i)-mrtab(i-1))*vrot(i)
        

        mvert=abs(zglass(nn))

        j=1
     
        do
           if(j.ge.nzbin-1) exit
           if(mztab(i-1,j).ge.mvert) exit
           j=j+1
        end do
        
        z=ztab(j)
        if(zglass(nn).lt.0) z=-1*z

        if(j.gt.nzbin-20) write(33,*) z,j,nzbin,mztab(i-1,j),mvert,zglass(nn)
        
        pos(1,nn)=r1*xglass(nn)
        pos(2,nn)=r1*yglass(nn)
        pos(3,nn)=z

        xc=xc+pos(1,nn)
        yc=yc+pos(2,nn)
        zc=zc+pos(3,nn)
        
        rmax=max(rmax,sqrt(pos(1,nn)**2+pos(2,nn)**2))
        
        vel(1,nn)=-vcirc*sin(atan2(pos(2,nn),pos(1,nn)))
        vel(2,nn)=vcirc*cos(atan2(pos(2,nn),pos(1,nn)))
        vel(3,nn)=0.0

        vxc=vxc+vel(1,nn)
        vyc=vyc+vel(2,nn)
        vzc=vzc+vel(3,nn)
        
        vxc2=vxc2+vel(1,nn)**2
        vyc2=vyc2+vel(2,nn)**2
        vzc2=vzc2+vel(3,nn)**2     
        
        id(nn)=nn
        
        ug(nn)=cs*cs/(gamma-1)

        jx=jx+(pos(2,nn)*vel(3,nn)-pos(3,nn)*vel(2,nn))
        jy=jy-(pos(1,nn)*vel(3,nn)-pos(3,nn)*vel(1,nn))
        jz=jz+(pos(1,nn)*vel(2,nn)-pos(2,nn)*vel(1,nn))
     end do
     deallocate(xglass,yglass,zglass)
  end if

  xc=xc/real(nn)
  yc=yc/real(nn)
  zc=zc/real(nn)

  vxc=vxc/real(nn)
  vyc=vyc/real(nn)
  vzc=vzc/real(nn)

  vxc2=vxc2/real(nn)
  vyc2=vyc2/real(nn)
  vzc2=vzc2/real(nn)
  
  jx=jx/real(nn)
  jy=jy/real(nn)
  jz=jz/real(nn)    

  write(*,*) 'Mass of gas within disc [1e10 solar masses] : ',real(nn)*mpart
  write(*,*) 'Particle mass [1e10 Msol]:',mpart
  write(*,*) 'Centre of mass [kpc]: ',xc,yc,zc
  write(*,*) 'Centre of mass velocity [km/s]: ',vxc,vyc,vzc
  write(*,*) 'Centre of mass velocity dispersion [km/s]: ',sqrt(vxc2),sqrt(vyc2),sqrt(vzc2)
  write(*,*) 'Angular momentum [kpc x km/s]: ',jx,jy,jz,sqrt(jx*jx+jy*jy+jz*jz)/sqrt(2.*ggdt*m200*r200)  
  write(*,*)

  ! Measure the radial profiles we've just created...
  if(allocated(r2)) deallocate(r2)
  if(allocated(indx)) deallocate(indx)  
  allocate(r2(nn),indx(nn))
  
  do i=1,nn
     r2(i)=0.0
     do j=1,2
        r2(i)=r2(i)+pos(j,i)**2
     end do
  end do
  
  call indexx(nn,r2,indx)
  
  nring=0
  rring=0.0
  mring=0.0
  m_enc=0.0

  xmin=1.e10; ymin=1.e10; zmin=1.e10
  xmax=-1.e10; ymax=-1.e10; zmax=-1.e10
  xc=0.0; yc=0.0; zc=0.0; zc2=0.0; 
  jx=0.0; jy=0.0; jz=0.0
  
  afac=(1/0.01)**(0.05)
  rmax=0.001
  r0j=0.0

  open(2,file='./disc_radial_profile.txt',status='unknown')
  write(2,*) m200,r200,mdisc,mstar
  do i=1,nn
     j=indx(i)
     rj=sqrt(r2(j))
     nring=nring+1
     rring=rring+rj
     
     xmin=min(xmin,pos(1,j))
     xmax=max(xmax,pos(1,j))
     ymin=min(ymin,pos(2,j))     
     ymax=max(ymax,pos(2,j))
     zmin=min(zmin,pos(3,j))
     zmax=max(zmax,pos(3,j))

     zc2=zc2+pos(3,j)**2

     jx=jx+(pos(2,j)*vel(3,j)-pos(3,j)*vel(2,j))
     jy=jy-(pos(1,j)*vel(3,j)-pos(3,j)*vel(1,j))
     jz=jz+(pos(1,j)*vel(2,j)-pos(2,j)*vel(1,j))
     
     mring=mring+mpart
     m_enc=m_enc+mpart
     
     if(rj/r1.gt.rmax) then
        rring=rring*(1.0/dble(nring))
        sigma=mring/(PI*(rj**2.0-r0j**2.0))
        zc2=sqrt(zc2/real(nring))
        jx=jx/real(nring)
        jy=jy/real(nring)
        jz=jz/real(nring)        
        write(2,'(2e18.8,i9,6e18.8)') rring,sigma,nring,mring,m_enc,zc2,jx,jy,jz
        jx=0.0; jy=0.0; jz=0.0; zc2=0.0
        r0j=rj
        nring=0
        rring=0.0
        
        mring=0.0
        rmax=rmax*afac
     end if
  end do
  close(2)

  ! Set particle species
  
  do i=1,6
     np(i)=0
     massarr(i)=0.0
     highword(i)=0.0
  end do
  
  np(1)=nn
  massarr(1)=mpart
  
  if(isbh.eqv..true.) then
     np(6)=1
     massarr(6)=mbh
     pos(1,nn+1)=0.0
     pos(2,nn+1)=0.0
     pos(3,nn+1)=0.0
     
     vel(1,nn+1)=0.0
     vel(2,nn+1)=0.0
     vel(3,nn+1)=0.0
  end if

  do i=1,6
     nall(i)=np(i)
     highword(i)=0!ishft(np(i),-16)
  end do
  
  ! Write out data

#ifdef HDF5  
  if(mode.eq.1.or.mode.eq.2) then
#endif
     call write_gadget_bin(outfile,np(1)+np(6),np(1),mode)
#ifdef HDF5  
  else
     call write_hdf5(outfile,np(1)+np(6),np(1),mode)     
  end if
#endif
  
  deallocate(pos,vel,id,ug)
  
end program main
    
    
real(kind=8) function dint(k,z)
  use structure
  use nrutils_modules
  implicit none
  
  real(kind=8) :: k,z
  
  dint = bessj0(k*radius)*dexp(-k*z)
  dint = dint/(1+(k*rdisc)**2)**1.5
  
end function dint

real(kind=8) function fint(x)
  use structure
  implicit none
      
  real(kind=8) :: x
  
  real(kind=8), external :: phiz
  
  fint = exp(-phiz(x)/cs**2)
  
end function fint

real(kind=8) function phiz(z)
  use constants
  use structure
  use test
  implicit none
  
  real(kind=8) :: z  ! Projected radius, height above midplane
  real(kind=8) :: fofc,phi,lphi
  
  real(kind=8) :: x,a,mtotal
  
  integer(kind=4) :: i,j
  
  i=1
  do
     if(rtab(i).ge.radius) exit
     i=i+1
  end do
  
  j=1
  do
     if(ztab(j).ge.z) exit
     j=j+1
  end do
  
  phiz=0
  
  if(j.eq.1) return  ! In this case, makes no difference what radial position is - we are at base
  ! of midplane and so potential difference wrt z is zero
  
  if(i.eq.1) then
     phi = pot(1,1)
     phiz = (z-ztab(j-1))/(ztab(j)-ztab(j-1))*pot(1,j)+&
          & (ztab(j)-z)/(ztab(j)-ztab(j-1))*pot(1,j-1)
     phiz= phiz-phi
     return
  end if
  
  phi=(radius-rtab(i-1))/(rtab(i)-rtab(i-1))*pot(i,1)+&
       & (rtab(i)-radius)/(rtab(i)-rtab(i-1))*pot(i-1,1)
  
  phiz=(radius-rtab(i-1))/(rtab(i)-rtab(i-1))*&
       &    ((z-ztab(j-1))/(ztab(j)-ztab(j-1))*pot(i,j)+&
       &    (ztab(j)-z)/(ztab(j)-ztab(j-1))*pot(i,j-1))+&
       &    (rtab(i)-radius)/(rtab(i)-rtab(i-1))*&
       &    ((z-ztab(j-1))/(ztab(j)-ztab(j-1))*pot(i-1,j)+&
       &    (ztab(j)-z)/(ztab(j)-ztab(j-1))*pot(i-1,j-1))
  
  phiz=phiz-phi
  
  return
  
end function phiz

subroutine derivs(x,y,dydx)
  use constants
  use structure
  implicit none
  
  real(kind=8), intent(in) :: x,y(2)
  real(kind=8), intent(out) :: dydx(2)
  
  real(kind=8), external :: phiz
  
  dydx(1)=y(2)
  dydx(2)=4.*PI*ggdt*rho0*dexp(-(phiz(x)+y(1))/cs**2)
end subroutine derivs


real(kind=8) function fintrp(x,xp,yp,np)
  use constants
  use structure
  implicit none
  
  real(kind=8), external :: phiz
  
  integer(kind=4) :: n,np,nhi,nlo
  real(kind=8), dimension(np) :: xp,yp
  real(kind=8) :: x,phigz,dx,x1,xn,h
  
  n=1
  
  do
     if(xp(n).ge.x) exit
     n=n+1
  end do
  
  if(xp(n).eq.x) then
     phigz=yp(n)
     fintrp = exp(-(phiz(x)+phigz)/cs**2)
     return
  end if
  
  dx=xp(n)-xp(n-1)
  
  phigz = yp(n-1)*(xp(n)-x)+yp(n)*(x-xp(n-1))
  phigz = phigz/dx
  
  fintrp = exp(-(phiz(x)+phigz)/cs**2)
  
  return
  
end function fintrp

