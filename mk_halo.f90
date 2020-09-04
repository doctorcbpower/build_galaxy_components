program mk_halo
  use nrutils_modules
  use gadget_header
  use structure
  use read_params
  use gas_halo
  use io
  use make_input  
  
  implicit none
  
  integer(kind=4) :: npart,nkeep,ngvir,nbh,nsample
  
  ! Input data
  integer(kind=4) :: nglass,nexp,ntot
  real(kind=4), allocatable :: xglass(:), yglass(:), zglass(:)

  real(kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=8) :: xc,yc,zc,vxc,vyc,vzc,vxc2,vyc2,vzc2,lbox
  real(kind=4), allocatable :: r2(:)
  integer(kind=4), allocatable :: indx(:)

  real(kind=8) :: jx,jy,jz,jtot
  ! Output data
  real(kind=4) :: mpart,logmp,msample,fmsample
  
  real(kind=8) :: dx,dy,dz,r,rvir,rcynd,phi

  character(kind=1,len=20) :: instring
  character(kind=1,len=132) :: infile,outfile,datafile,glassfile,paramfile
  logical :: fexist

  integer(kind=4) :: i,ii,j,k,n

  real(kind=8) :: vshell
  real(kind=8) :: utemp,rhobar,vesc
  real(kind=8) :: lx,lxmin,lxmax,dlx

  integer(kind=4), parameter :: ntab=1000
  real(kind=8), parameter :: fmin=1.0e-5
  real(kind=8), dimension(ntab) :: rtab,mtab,utab
  real(kind=8) :: afac,mrmin,mrvir,mrmax,mtotal,rmin,rmax,mgas

  real(kind=8) :: rnew

  integer(kind=4) :: nshell
  real(kind=8) :: rj,r0j
  real(kind=8) :: rshell,mshell,m_enc,menclo,rho,rsc,router

  real(kind=8) :: a
  logical :: isbh

  logical :: ispoisson

  integer(kind=4) :: mode
  isbh=.false.

  nexp=3

#if defined(HERNQUIST_HALO) && defined(NFW_HALO)
  write(*,*) 'Error: multiple DM halos!'
#endif      

  if(command_argument_count().eq.0) stop 'Usage: mk_halo.exe ./paramfile'
  
  call get_command_argument(1,paramfile)

  inquire(file=paramfile,exist=fexist)

  if(fexist.eqv..false.) stop 'Error: parameter file does not exist'
  call set_parameters(paramfile,npart,outfile,glassfile,ispoisson,1,mode)
  if(mbh.gt.0.0) isbh=.true.

  write(*,*)
  write(*,*) 'Defining virial overdensity (wrt critical density) as: ',deltavir
  write(*,*)

  if(v200.gt.0.0) then
     mvir=getmvir(v200,deltavir,rhocrit0)
     if(mvir.ne.m200) m200=mvir
  end if
  
  r200 = getrvir(m200,deltavir,rhocrit0)
  rs=r200/c200                                  ! NFW Scale radius
  router=1.1*r200
  
  if(isbh.eqv..true.) write(*,*) 'Assuming a central black hole...'

#ifdef HERNQUIST_HALO
  write(*,*) 'Assuming a Hernquist halo...'
  write(*,*) 'Halo properties:'
  write(*,*) 'Mvir: ',m200
  write(*,*) 'Rvir: ',r200      
#ifdef GalIC
  write(*,*) 'Mtot: ',m200*(1+(1/c200)*sqrt(2*(dlog(1+c200)-c200/(1+c200))))
#else
  write(*,*) 'Mtot: ',m200
#endif
  write(*,*) 'Scale radius:',rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
  write(*,*)
#endif
  
#ifdef NFW_HALO
  write(*,*) 'Assuming a NFW halo...'
  write(*,*) 'Halo properties:'
  write(*,*) 'Mvir: ',m200
  write(*,*) 'Rvir: ',r200      
  write(*,*) 'Scale radius:',rs
  write(*,*)
#endif            
  
!!! Now set up the enclosed mass profile 

  write(*,*) 'Assuming a baryon fraction :',fbaryon

  xmin=1.0e-6             ! These distances are in units of r200
  xmax=1.1*router/r200

  call menc_gas(xmin,xmax,mgas,fbaryon,ntab,rtab,mtab,utab)

  ! Want to define the particle mass to ensure nhalo particles inside r200
  logmp=log10(mgas/real(npart))

#ifdef HERNQUIST_HALO
  a=rs*sqrt(2.*(dlog(1+c200)-c200/(1+c200)))
  
  i=1
  
  do 
     if(rtab(i).gt.router) exit
     i=i+1
  end do

  msample=mgas*(mtab(i-1)*(rtab(i)-router)/(rtab(i)-rtab(i-1))+&
       & mtab(i)*(router-rtab(i-1))/(rtab(i)-rtab(i-1)))
#endif
  
#ifdef NFW_HALO  
  i=1
  
  do 
     if(rtab(i).gt.router) exit
     i=i+1
  end do

  msample=mgas*(mtab(i-1)*(rtab(i)-router)/(rtab(i)-rtab(i-1))+&
          & mtab(i)*(router-rtab(i-1))/(rtab(i)-rtab(i-1)))
#endif
  
  nsample=1+floor(msample*(10**(-logmp)))
  fmsample=msample/mgas
  
  if(ispoisson.eqv..false.) then
     inquire(file=glassfile,exist=fexist)
     if(fexist.eqv..false.) stop 'Could not find input glass file...'
     call get_nglass(glassfile,nglass)

     ! nsample is the number of particles within a sphere - want to
     ! estimate how many particles in the corresponding cube with the
     ! same number density is required

     nexp=1+floor((6*real(nsample)/pi/real(nglass))**(1./3.))

     ntot=nexp*nexp*nexp*nglass

     allocate(xglass(ntot))
     allocate(yglass(ntot))
     allocate(zglass(ntot))
     if(allocated(xglass).eqv..false. .or.&
          & allocated(yglass).eqv..false. .or.&
          & allocated(zglass).eqv..false.)&
          &  stop 'Could not allocate memory!'     
     call get_glass_particles(glassfile,0,nexp,nglass,ntot,xglass,yglass,zglass,0)
  else
     write(*,*) 'Sampling ',nsample,' particles from a Poisson distribution'     
  end if
  
  nglass=nsample

  allocate(pos(3,nglass+1))
  if(allocated(pos).eqv..false.) stop 'Error allocating memory (pos)!'
  
  allocate(vel(3,nglass+1))
  if(allocated(vel).eqv..false.) stop 'Error allocating memory (vel)!'

  allocate(id(nglass+1))
  if(allocated(id).eqv..false.) stop 'Error allocating memory (id)!'
  
  allocate(ug(nglass))
  if(allocated(ug).eqv..false.) stop 'Error allocating memory (ug)!'

  mtotal=0.0

  nkeep = 0
  ngvir = 0
  vshell = 0.0
  
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

  if(ispoisson) iseed=-18827921
  
  do 
     if(nkeep.eq.nglass) exit
     
     nkeep=nkeep+1
     
     if(ispoisson.eqv..false.) then
        dx=xglass(nkeep)-0.5
        dy=yglass(nkeep)-0.5
        dz=zglass(nkeep)-0.5
        
        r=sqrt(dx*dx+dy*dy+dz*dz)
        
        m_enc=real(nkeep)/real(nglass)*fmsample
     else
#ifdef NFW_HALO     
        m_enc = ran3(iseed)*fmsample
#endif
#ifdef HERNQUIST_HALO     
        m_enc = ran3(iseed)*fmsample
#endif     
     end if
     
     i=1
     
     do 
        if(mtab(i).gt.m_enc) exit
        i=i+1
     end do
     
     if(i.eq.1) stop 'Need to extend lower limit'

     rnew=rtab(i-1)*(mtab(i)-m_enc)/(mtab(i)-mtab(i-1))+&
          & rtab(i)*(m_enc-mtab(i-1))/(mtab(i)-mtab(i-1))

     if(ispoisson) then
        phi=2*PI*ran3(iseed)
        
        pos(1,nkeep)=rnew*(2.0*ran3(iseed)-1.0)
        rcynd=dsqrt(rnew*rnew-pos(1,nkeep)*pos(1,nkeep))
        pos(2,nkeep)=rcynd*cos(phi)
        pos(3,nkeep)=rcynd*sin(phi)
     else
        pos(1,nkeep)=rnew*dx/r
        pos(2,nkeep)=rnew*dy/r
        pos(3,nkeep)=rnew*dz/r
     end if
     
     xc=xc+pos(1,nkeep)
     yc=yc+pos(2,nkeep)
     zc=zc+pos(3,nkeep)
     
     rmax=max(rmax,sqrt(pos(1,nkeep)**2+pos(2,nkeep)**2+pos(3,nkeep)**2))
     
     if(sqrt(pos(1,nkeep)**2+pos(2,nkeep)**2+pos(3,nkeep)**2).le.r200)&
          & ngvir=ngvir+1

     ! Include the radial component...
     
     vel(1,nkeep) = v_radial * (pos(1,nkeep)/rnew) * vunit
     vel(2,nkeep) = v_radial * (pos(2,nkeep)/rnew) * vunit
     vel(3,nkeep) = v_radial * (pos(3,nkeep)/rnew) * vunit

     ! ... and angular component if present

     vel(1,nkeep) = vel(1,nkeep) + get_vel_comps(pos(1:3,nkeep),1)
     vel(2,nkeep) = vel(2,nkeep) + get_vel_comps(pos(1:3,nkeep),2)
     vel(3,nkeep) = vel(3,nkeep) + get_vel_comps(pos(1:3,nkeep),3)

     vxc=vxc+vel(1,nkeep)
     vyc=vyc+vel(2,nkeep)
     vzc=vzc+vel(3,nkeep)

     vxc2=vxc2+vel(1,nkeep)**2
     vyc2=vyc2+vel(2,nkeep)**2
     vzc2=vzc2+vel(3,nkeep)**2     
     
     id(nkeep)=nkeep
     
     ug(nkeep) = utab(i-1)*(rtab(i)-rnew)/(rtab(i)-rtab(i-1))+&
          & utab(i)*(rnew-rtab(i-1))/(rtab(i)-rtab(i-1))
     ug(nkeep) = ug(nkeep)*vunit**2.

     jx=jx+(pos(2,nkeep)*vel(3,nkeep)-pos(3,nkeep)*vel(2,nkeep))
     jy=jy-(pos(1,nkeep)*vel(3,nkeep)-pos(3,nkeep)*vel(1,nkeep))
     jz=jz+(pos(1,nkeep)*vel(2,nkeep)-pos(2,nkeep)*vel(1,nkeep))
     
     mtotal=mtotal+10**logmp
  end do

  mpart = 10**logmp
  
  xc=xc/real(nkeep)
  yc=yc/real(nkeep)
  zc=zc/real(nkeep)

  vxc=vxc/real(nkeep)
  vyc=vyc/real(nkeep)
  vzc=vzc/real(nkeep)

  vxc2=vxc2/real(nkeep)
  vyc2=vyc2/real(nkeep)
  vzc2=vzc2/real(nkeep)
  
  jx=jx/real(nkeep)
  jy=jy/real(nkeep)
  jz=jz/real(nkeep)    
  
  write(*,*) 'Mass of gas within rvir [1e10 solar masses] : ',real(ngvir)*mpart
  write(*,*) 'Mass of gas within router [1e10 solar masses] : ',real(nkeep)*mpart
  write(*,*) 'Baryon fraction : ',real(ngvir)*mpart/m200
  write(*,*) 'Particle mass [1e10 Msol]:',mpart
  write(*,*) 'Centre of mass [kpc]: ',xc,yc,zc
  write(*,*) 'Centre of mass velocity [km/s]: ',vxc,vyc,vzc
  write(*,*) 'Centre of mass velocity dispersion [km/s]: ',sqrt(vxc2),sqrt(vyc2),sqrt(vzc2)
  write(*,*) 'Angular momentum [kpc x km/s]: ',jx,jy,jz,sqrt(jx*jx+jy*jy+jz*jz)/sqrt(2.*ggdt*m200*r200)  

  ! Measure the radial profiles we've just created...
  if(allocated(r2)) deallocate(r2)
  if(allocated(indx)) deallocate(indx)  
  allocate(r2(nkeep),indx(nkeep))

  do i=1,nkeep
     r2(i)=0.0
     do j=1,3
        r2(i)=r2(i)+pos(j,i)**2
     end do
  end do
  
  call indexx(nkeep,r2,indx)
  
  nshell=0
  rshell=0.0
  mshell=0.0
  m_enc=0.0

  xmin=1.e10; ymin=1.e10; zmin=1.e10
  xmax=-1.e10; ymax=-1.e10; zmax=-1.e10
  xc=0.0; yc=0.0; zc=0.0
  
  afac=(1/0.001)**(0.0025)
  rmax=0.001
  r0j=0.0

  open(2,file='./halo_radial_profile.txt',status='unknown')
  write(2,*) m200,r200,mgas,fbaryon
  do i=1,nkeep
     j=indx(i)
     rj=sqrt(r2(j))
     nshell=nshell+1
     rshell=rshell+rj

     xmin=min(xmin,pos(1,j))
     xmax=max(xmax,pos(1,j))
     ymin=min(ymin,pos(2,j))     
     ymax=max(ymax,pos(2,j))
     zmin=min(zmin,pos(3,j))
     zmax=max(zmax,pos(3,j))

     if(i.lt.nkeep/2) then
        xc=xc+pos(1,j)
        yc=yc+pos(2,j)
        zc=zc+pos(3,j)
     end if
     
     mshell=mshell+mpart
     m_enc=m_enc+mpart
     
     if(rj/r200.gt.rmax) then
        rshell=rshell*(1.0/dble(nshell))
        rho=mshell/(4.*PI*(rj**3.0-r0j**3.0)/3.)
        write(2,'(2e18.8,i9,2e18.8)') rshell,rho,nshell,mshell,m_enc
        
        r0j=rj
        nshell=0
        rshell=0.0
        
        mshell=0.0
        rmax=rmax*afac
     end if
  end do
  close(2)

  xc=2*xc/real(nkeep); yc=2*yc/real(nkeep); zc=2*zc/real(nkeep)

  lbox=1.+floor(max(xmax-xmin,ymax-ymin,zmax-zmin))

  do i=1,nkeep
     pos(1,i)=pos(1,i)-xc
     pos(2,i)=pos(2,i)-yc
     pos(3,i)=pos(3,i)-zc
  end do
  
  ! Now write the data to file...
  
  ! Set particle species
  
  do i=1,6
     np(i)=0
     massarr(i)=0.0
     highword(i)=0.0
  end do
  
  np(1)=nkeep
  massarr(1)=mpart
  
  if(isbh.eqv..true.) then
     np(6)=1
     massarr(6)=mbh
     pos(1,nkeep+1)=0.0
     pos(2,nkeep+1)=0.0
     pos(3,nkeep+1)=0.0
     
     vel(1,nkeep+1)=0.0
     vel(2,nkeep+1)=0.0
     vel(3,nkeep+1)=0.0
  end if

  BoxSize=lbox
  
  do i=1,6
     nall(i)=np(i)
!     highword(i)=ishft(np(i),-16)
  end do

#ifdef HDF5    
  if(mode.eq.1.or.mode.eq.2) then
#endif
     call write_gadget_bin(outfile,nall(1)+nall(6),nall(1),mode)
#ifdef HDF5  
  else
     call write_hdf5(outfile,nall(1)+nall(6),nall(1),mode)
  end if
#endif  
  deallocate(pos,vel,id,ug)
  
end program mk_halo



  
