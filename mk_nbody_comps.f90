module params
  implicit none
  real(kind=8) :: xx,a
end module params

module prof_params
  implicit none

  real(kind=8) :: rho_0    ! Characteristic density of halo
  real(kind=8) :: pot_0    ! Characteristic density of halo
  real(kind=8) :: pot_in,pot_out
  real(kind=8) :: rvir     ! Scale radius
  real(kind=8) :: rscale   ! Scale radius

  integer(kind=4) :: icomp  ! Component index

  integer(kind=4), parameter :: ntab=1000
  integer(kind=4) :: nkeep
  real(kind=8), dimension(ntab) :: rtab,ptab,mtab,dptab,drptab,dddtab,etab
  real(kind=8), dimension(5,ntab) ::  dtab,drtab,ddrtab,drhopsitab,ftab,rhotab,menc_tab
  real(kind=8) :: b_t,rstar,rho_star,mu  ! Parameters needed to taper mass profiles
  real(kind=8) :: r_smooth, z_smooth
end module prof_params

program mk_halo
  
  use params
  use prof_params
  use nrutils_modules
  use constants
  use structure
  use read_params
  use io

  implicit none

  real(kind=8) :: lemin,lemax,dle,le
  real(kind=8) :: lrmin,lrmax,dlr,lr

  real(kind=8) :: rmax
  
  real(kind=8) :: mtot,fsum,fmax,eran,fran,fsamp,msum,mtotal,sigma
  real(kind=8) :: mtot_in_comp,mtot_in_m200
  
  real(kind=8), external :: fabel,fpot

  integer(kind=4) :: i,j,k,l,m,n,nmin,nmax

  real(kind=8) :: xmin,xmid,xmax
  real(kind=8) :: rcynd,phi,theta,mofr,vmag,vesc,mpart
  real(kind=8) :: rav,r0,v1,v2,v3,vrav,x1,x2,x3
  real(kind=8) :: vxav,vyav,vzav,vx2av,vy2av,vz2av,vr2av
  real(kind=8) :: xav,yav,zav
  real(kind=8), parameter :: acc=1d-6

  integer(kind=4) :: ncount,npart,npart_comp_tot,nvir

  ! Process command line arguments
  character(len=132), dimension(20) :: instring ! command line arguments

  integer(kind=4) :: imvir,iprof,iconc,inpart,imax,ifile,iran,isigma,isoft

  integer(kind=4) :: nit
  integer(kind=4), parameter :: nitmax=1000

  real(kind=4) :: xoff,yoff,zoff,vxoff,vyoff,vzoff
  real(kind=8) :: ell,rp,ra,vp,vcirc,softening

  real(kind=8) :: dlrhodlrl,dlrhodlru
  real(kind=8), dimension(ntab) :: logr,logrho,d2lrhodlr2
  
  real(kind=8), dimension(5) :: rho_comp,mass_comp,f_comp
  integer(kind=4) ::n_in_comp,ncmp
  real(kind=8) :: menc_r

  real(kind=8), external :: dmass,fmass

!  real(kind=8), external :: inner_int

  character(kind=1,len=132) :: paramfile,glassfile,outfile,tabfile
  logical:: fexist
  logical:: ispoisson
  logical:: istabulated  
  integer(kind=4) :: mode           ! mode of output
  logical:: isbh

  real(kind=8), parameter :: frcut=1.05
  real(kind=8) :: mrcut

  isbh=.false.
  
  if(command_argument_count().eq.0) stop 'Usage: mk_nbody_cpmp.exe <parameter_file>'

  call get_command_argument(1,paramfile)

  inquire(file=paramfile,exist=fexist)

  if(fexist.eqv..false.) stop 'Error: parameter file does not exist'
  call set_parameters(paramfile,npart,outfile,glassfile,istabulated,tabfile,ispoisson,2,mode)

  if(mbh.gt.0.0) isbh=.true.

  if(isbh.eqv..true.) write(*,*) 'Assuming a central black hole...'

  do i=1,ncomp
     if(raniso_comp(i).gt.0.0) write(*,'(a,i2)') 'Assuming an anistropic velocity distribution for component',i
  end do
  
  ! Compute properties of the halo
  
  if(v200.gt.0.0) then
    mvir=getmvir(v200,deltavir,rhocrit0)
    if(mvir.ne.m200) m200=mvir
  end if
 
  r200=getrvir(m200,deltavir,rhocrit0)
  rs=r200/c200                                  ! NFW Scale radius

  write(*,*) 'Virial radius [kpc/h]:',r200
  write(*,*) 'Scale radius [kpc/h]:',rs

  lrmin=dlog10(rs)-6.
  lrmax=dlog10(rs)+6.
  
  if(istaper) call taper_mass_profile(mtot,10**lrmax)

  do i=1,ncomp
     prof_type=1
     if(form_comp(i).eq.'hernquist') prof_type=2
     if(form_comp(i).eq.'burkert') prof_type=3
     if(form_comp(i).eq.'general') then
        prof_type=4
        alpha=1.
        beta=4.
        gamma_slope=2.
        epsilon=(-gamma_slope-beta*c200**alpha)/(1.+c200**alpha)+1.0/fdecay
     end if
     mass_comp(i)=fmass_comp(i)*m200
     rscale_comp(i)=rscale_comp(i)*rs
     rho_comp(i)=get_rho0(mass_comp(i),r200,rscale_comp(i))
  end do

  ! Now tabulate the mass profile -- allows for a multi-component system

  lrmin=lrmin*ln10   ! Convert to natural logarithms
  lrmax=lrmax*ln10
  
  dlr=(lrmax-lrmin)/real(ntab-1)
  lr=lrmin

  n=0

  do
     if((lr-lrmax).gt.1.0d-6) exit
     n=n+1
     rtab(n)=dexp(lr)
     do icomp=1,ncomp
        prof_type=1
        if(form_comp(icomp).eq.'hernquist') prof_type=2
        if(form_comp(icomp).eq.'burkert') prof_type=3
        if(form_comp(icomp).eq.'general') prof_type=4                
        if(istaper) then
           rhotab(icomp,n)=rho_trunc(rtab(n),rscale_comp(icomp),rho_comp(icomp),&
                & r200*rtrunc_comp(icomp),rstar_comp(icomp),rhostar_comp(icomp))
        else
           mhalo = mass_comp(icomp) 
           rs=rscale_comp(icomp)
           rhotab(icomp,n)=rhodm(rtab(n))
        end if
        menc_tab(icomp,n)=0.0
     end do
     mtab(n)=0.0
     lr=lr+dlr
  end do

  if(isbh.eqv..true.) mtab(1)=mbh  

  do icomp=1,ncomp
     prof_type=1
     if(form_comp(icomp).eq.'hernquist') prof_type=2
     if(form_comp(icomp).eq.'burkert') prof_type=3
     if(form_comp(icomp).eq.'general') prof_type=4                     
     rs = rscale_comp(icomp)
     mhalo = mass_comp(icomp) 
     menc_tab(icomp,1)=mhenc(rtab(1))
     mrcut=0.0
     do n=1,ntab
        call qromo(dmass,lrmin,dlog(rtab(n)),menc_r,midpnt)
        if(istaper.and.rtab(n).ge.frcut*r200*rtrunc_comp(icomp)) then
           if(mrcut.eq.0.0) mrcut=menc_r
           menc_r=mrcut
        end if
        menc_tab(icomp,n)=menc_tab(icomp,1)+menc_r
        mtab(n)=mtab(n)+menc_tab(icomp,n)
     end do
  end do

  i=0

  do
     if(i.gt.ntab) stop
     i=i+1
     if(rtab(i).gt.r200) exit
  end do

  mtot=mtab(i-1)*(rtab(i)-r200)/(rtab(i)-rtab(i-1))+&
       & mtab(i)*(r200-rtab(i-1))/(rtab(i)-rtab(i-1))

  rscale = rscale_comp(1)
  
  write(*,*) 'Total mass [1e10 Msol/h]:',mtot 
  
  write(*,*) 'Characteristic potential :',log10(ggdt*mtot/rscale)

  write(*,*) 'Characteristic density :',log10(rho_comp(1))  
  
  write(*,*) 'Circular speed at virial radius [km/s] :',sqrt(ggdt*m200/r200)

  write(*,*)

  if(form_comp(1).eq.'nfw'.or.form_comp(1).eq.'burkert'.or.form_comp(1).eq.'general') then
     mtot = mass_comp(1)
  else
     mtot = mtab(ntab)
  end if
  
  pot_0=0.0

  do i=ntab,1,-1
     call qromo(fpot,dlog(rtab(i)),lrmax,ptab(i),midpnt)
     ptab(i)=-ggdt*ptab(i)
     if(i.eq.ntab) pot_0=ptab(i)
     ptab(i)=-(ptab(i)-pot_0)*rscale/ggdt/mtot
     print *, rtab(i),mtab(i),ptab(i)
  end do

  write(*,*) 'Tabulated gravitational potential'

  ! Now tabulate the density and potential, for inversion

  
  if((istaper.eqv..false.).and.(fdecay.eq.0.0)) then
     do icomp=1,ncomp
        do n=1,ntab
           dtab(icomp,n)=dlog(rhotab(icomp,n))
        end do
     end do
    
     ! Compute the 1st derivative of density with respect to radius
     ! Note that this is dln(rho)/dln(r)
     do n=3,ntab-2
        do icomp=1,ncomp
           drtab(icomp,n)=(-dtab(icomp,n+2)+8*dtab(icomp,n+1)-&
                & 8*dtab(icomp,n-1)+dtab(icomp,n-2))/12./dlr
        end do
     end do
    
     do icomp=1,ncomp
        drtab(icomp,1)=drtab(icomp,3)-2*(drtab(icomp,4)-drtab(icomp,3))
        drtab(icomp,2)=drtab(icomp,3)-(drtab(icomp,4)-drtab(icomp,3))
        drtab(icomp,ntab-1)=drtab(icomp,ntab-2)+(drtab(icomp,ntab-2)-drtab(icomp,ntab-3))
        drtab(icomp,ntab)=drtab(icomp,ntab-2)+2.*(drtab(icomp,ntab-2)-drtab(icomp,ntab-3))     
     end do
  else
     ! Compute the 1st derivative of density with respect to radius
     ! Note that this is drho/dln(r)
     do n=3,ntab-2
        do icomp=1,ncomp
           drtab(icomp,n)=(-rhotab(icomp,n+2)+8*rhotab(icomp,n+1)-&
                & 8*rhotab(icomp,n-1)+rhotab(icomp,n-2))/12./dlr
        end do
     end do
     
     do icomp=1,ncomp
        drtab(icomp,1)=drtab(icomp,3)-2*(drtab(icomp,4)-drtab(icomp,3))
        drtab(icomp,2)=drtab(icomp,3)-(drtab(icomp,4)-drtab(icomp,3))
        drtab(icomp,ntab-1)=drtab(icomp,ntab-2)+(drtab(icomp,ntab-2)-drtab(icomp,ntab-3))
        drtab(icomp,ntab)=drtab(icomp,ntab-2)+2.*(drtab(icomp,ntab-2)-drtab(icomp,ntab-3))     
     end do
  end if

  ! Compute the 2nd derivative of density with respect to radius
  do n=3,ntab-2
     do icomp=1,ncomp
        ddrtab(icomp,n)=(-drtab(icomp,n+2)+8*drtab(icomp,n+1)-&
             & 8*drtab(icomp,n-1)+drtab(icomp,n-2))/12./dlr
     end do
  end do
  
  do icomp=1,ncomp
     ddrtab(icomp,1)=ddrtab(icomp,3)-2.*(ddrtab(icomp,4)-ddrtab(icomp,3))
     ddrtab(icomp,2)=ddrtab(icomp,3)-(ddrtab(icomp,4)-ddrtab(icomp,3))
     ddrtab(icomp,ntab-1)=ddrtab(icomp,ntab-2)+(ddrtab(icomp,ntab-2)-ddrtab(icomp,ntab-3))
     ddrtab(icomp,ntab)=ddrtab(icomp,ntab-2)+2.*(ddrtab(icomp,ntab-2)-ddrtab(icomp,ntab-3))
  end do
     
  ! Compute the 1st derivative of relative potential with respect to radius -
  ! this is simply the inverse of the acceleration
  do n=1,ntab
     dptab(n)=rtab(n)*rtab(n)/ggdt/mtab(n)
  end do

  do n=3,ntab-2
     drptab(n)=(-dlog(mtab(n+2))+8*dlog(mtab(n+1))&
          & -8*dlog(mtab(n-1))+dlog(mtab(n-2)))/12./dlr
  end do
  
  drptab(1)=drptab(3)-2.*(dptab(4)-dptab(3))
  drptab(2)=drptab(3)-(dptab(4)-dptab(3))     
  drptab(ntab-1)=drptab(ntab-2)+(dptab(ntab-2)-dptab(ntab-3))
  drptab(ntab)=drptab(ntab-2)+2.*(dptab(ntab-2)-dptab(ntab-3))

  do n=1,ntab
     drptab(n)=(2.-drptab(n))*rtab(n)/ggdt/mtab(n)
  end do

  if(istaper.or.(fdecay.gt.0.0)) then
     do icomp=1,ncomp
        do n=1,ntab
           drhopsitab(icomp,n)=(1.0/rtab(n)**2)*&
                & (ddrtab(icomp,n)-drtab(icomp,n))*&
                & dptab(n)*dptab(n)
           drhopsitab(icomp,n)=drhopsitab(icomp,n)+(1.0/rtab(n))*&
                & drtab(icomp,n)*dptab(n)*drptab(n)
           drhopsitab(icomp,n)=drhopsitab(icomp,n)*(ggdt*ggdt*mtot)/(rs*rs)
        end do
     end do
  else
     do icomp=1,ncomp
        do n=1,ntab
           drhopsitab(icomp,n)=(dexp(dtab(icomp,n))/rtab(n)**2)*&
                & (ddrtab(icomp,n)+drtab(icomp,n)**2-drtab(icomp,n))*&
                & dptab(n)*dptab(n)
           drhopsitab(icomp,n)=drhopsitab(icomp,n)+(dexp(dtab(icomp,n))/rtab(n))*&
                & drtab(icomp,n)*dptab(n)*drptab(n)
           drhopsitab(icomp,n)=drhopsitab(icomp,n)*(ggdt*ggdt*mtot)/(rs*rs)
        end do
     end do
  end if
  
  i=ntab
  do
     if(mtab(i).lt.0.95*mtot) exit
     i=i-1
  end do

  write(*,*) 'Computing distribution function...'

  nkeep=ntab-1
  
  lemin=dlog10(ptab(ntab-1))
  lemax=dlog10(ptab(1))

!  if(lemin.lt.dlog10(ptab(nkeep))) lemin=dlog10(ptab(nkeep))+0.1
  
  dle=(lemax-lemin)/real(ntab-1)
  le=lemin
  
  n=0

  open(33,file='distribution_function.txt',status='unknown')
  write(33,*) '# Binding Energy E [g*mtot/rscale]  F(E)'
  do 
     if(le.gt.lemax) exit
     
     xx=10**le
     n=n+1
     etab(n)=xx

     do icomp=1,ncomp
        call qromo(fabel,1.d-20,xx,fsum,midpnt)
        ftab(icomp,n)=fsum/sqrt(8.)/pi/pi
        if(ftab(icomp,n).le.1.e-12) ftab(icomp,n)=0.0
     end do
     
     write(33,'(10e18.8)') xx,(rscale_comp(1)**3*ftab(icomp,n),icomp=1,ncomp)
     le=le+dle
  end do
  close(33)

  nkeep=n

  write(*,*) 'Sampling the density profile...'

  ncount=0; nvir=0

  xav=0.0; yav=0.0; zav=0.0; vxav=0.0; vyav=0.0; vzav=0.0; vx2av=0.0; vy2av=0.0; vz2av=0.0; vrav=0.0; vr2av=0.0;

  npart=0
  
  do icomp=1,ncomp
     mpart=mass_comp(icomp)/real(npart_comp(icomp))
     mtot=menc_tab(icomp,ntab)
     npart=npart+1+floor(mtot/mpart)
  end do

  write(*,*) 'Allocating memory for ',npart,' particles...'
  
  allocate(pos(3,npart))
  if(allocated(pos).eqv..false.) stop 'Error allocating memory (pos)!'
  
  allocate(vel(3,npart))
  if(allocated(vel).eqv..false.) stop 'Error allocating memory (vel)!'

  allocate(id(npart))
  if(allocated(vel).eqv..false.) stop 'Error allocating memory (vel)!'
  
  write(*,*) 'Total mass enclosed [1e10 Msol]:',mtab(ntab)

  mtotal=0.0

  do i=1,ntype
     np(i)=0
     massarr(i)=0.0
  end do

  iseed=19721
  
  do icomp=1,ncomp
     mpart=mass_comp(icomp)/real(npart_comp(icomp))
     mtot=menc_tab(icomp,ntab)
!     j=1
!     do
!        if(rtab(j).gt.10*r200) exit
!        j=j+1
!     end do
!    mtot=menc_tab(icomp,j)
     rs=rscale_comp(icomp)
     npart_comp_tot = 1+floor(mtot/mpart)
     mtot_in_comp=0.0
     mtot_in_m200=0.0
     n_in_comp=0
     
     do
        if(n_in_comp.ge.npart_comp_tot) exit
        ! Choose a random mass
999     mofr=ran3(iseed)
        
        nmin=1
        nmax=ntab

        if((menc_tab(icomp,nmin)/mtot-mofr)*(menc_tab(icomp,nmax)/mtot-mofr).gt.0.0) &
             & stop '(ptab(1)-y)*(ptab(nmax)-y).gt.0.0'
        
        do
           n=(nmin+nmax)/2
           if((menc_tab(icomp,n)/mtot-mofr)*(menc_tab(icomp,nmax)/mtot-mofr).ge.0.0) then
              nmax=n
           else
              nmin=n
           end if
           
           if(abs(nmax-nmin).le.1) exit
        end do
        
       radius = rtab(nmin)*(menc_tab(icomp,nmax)-mofr*mtot)/(menc_tab(icomp,nmax)-menc_tab(icomp,nmin)) + &
           & rtab(nmax)*(mofr*mtot-menc_tab(icomp,nmin))/(menc_tab(icomp,nmax)-menc_tab(icomp,nmin))

       ! We know what the radius is -- now we need to estimate the relative potential...
       i=1
       do
          if(rtab(i).gt.radius) exit
          i=i+1
       end do
       
       phi=ptab(i+1)*(radius-rtab(i))/(rtab(i+1)-rtab(i))+ptab(i)*(rtab(i+1)-radius)/(rtab(i+1)-rtab(i))
       
       vesc=dsqrt(2.d0*phi)

       nit=0

       do
          if(nit.gt.nitmax) goto 999          
          
          nit=nit+1
          
          v1=(2.*ran3(iseed)-1.)*vesc
          v2=(2.*ran3(iseed)-1.)*vesc
          v3=(2.*ran3(iseed)-1.)*vesc

          eran=phi-0.5*(v1*v1+v2*v2+v3*v3)
          
          if(eran.le.0.0) cycle

          nmin=1
          nmax=nkeep

          if((etab(1)-eran)*(etab(nmax)-eran).gt.0.0) cycle
          
          do
             n=(nmin+nmax)/2
             if((etab(n)-eran)*(etab(nmax)-eran).ge.0.0) then
                nmax=n
             else
                nmin=n
             end if
             
             if(abs(nmax-nmin).le.1) exit
          end do
                             
          fsamp=ftab(icomp,nmax)*(eran-etab(nmin))/(etab(nmax)-etab(nmin))+&
               & ftab(icomp,nmin)*(etab(nmax)-eran)/(etab(nmax)-etab(nmin))

          ! We also need to estimate where the maximum of the distribution function lies...
          
          nmin=1
          nmax=nkeep          
          
          if((etab(1)-phi)*(etab(nmax)-phi).gt.0.0) cycle

          do
             n=(nmin+nmax)/2
             if((etab(n)-phi)*(etab(nmax)-phi).ge.0.0) then
                nmax=n
             else
                nmin=n
             end if
             
             if(abs(nmax-nmin).le.1) exit
          end do
          
          fmax=ftab(icomp,nmax)*(phi-etab(nmin))/(etab(nmax)-etab(nmin))+&
               & ftab(icomp,nmin)*(etab(nmax)-phi)/(etab(nmax)-etab(nmin))
          
          fran=fmax*ran3(iseed)
          
          if(fsamp.gt.fran) exit
          
       end do
       
2102   ncount=ncount+1
       
       ! Now assign the selected radii cartesian coordinates.
       ! These will be uniformly distributed on the unit sphere.

       phi=2*pi*ran3(iseed)
       
       do
          x1=(2.0*ran3(iseed)-1.0)
          rcynd=sqrt(1.0-x1*x1)
          x2=rcynd*cos(phi)
          x3=rcynd*sin(phi)
          if(abs(1.0-(x1*x1+x2*x2+x3*x3)).lt.acc) exit
       end do
       
       pos(1,ncount)=radius*x1
       pos(2,ncount)=radius*x2
       pos(3,ncount)=radius*x3       

       vel(1,ncount)=v1*sqrt(ggdt*mass_comp(1)/rscale_comp(1))
       vel(2,ncount)=v2*sqrt(ggdt*mass_comp(1)/rscale_comp(1))
       vel(3,ncount)=v3*sqrt(ggdt*mass_comp(1)/rscale_comp(1))       
       
       id(ncount)=ncount

       mtot_in_comp=mtot_in_comp+mpart
       
       if(radius.le.r200) then
          mtot_in_m200=mtot_in_m200+mpart
          nvir=nvir+1
          
          vxav=vxav+vel(1,ncount)
          vyav=vyav+vel(2,ncount)
          vzav=vzav+vel(3,ncount)
          
          vx2av=vx2av+vel(1,ncount)**2
          vy2av=vy2av+vel(2,ncount)**2
          vz2av=vz2av+vel(3,ncount)**2
          
          vrav=vrav+(vel(1,ncount)*pos(1,ncount)+vel(2,ncount)*pos(2,ncount)+vel(3,ncount)*pos(3,ncount))/radius
          vr2av=vr2av+((vel(1,ncount)*pos(1,ncount)+vel(2,ncount)*pos(2,ncount)+vel(3,ncount)*pos(3,ncount))/radius)**2
          
          xav=xav+pos(1,ncount)
          yav=yav+pos(2,ncount)
          zav=zav+pos(3,ncount)
       end if
       
       n_in_comp=n_in_comp+1
       mtotal=mtotal+mpart
    end do
    np(icomp+1)=n_in_comp
    massarr(icomp+1)=mpart
 end do

 xav=xav/real(nvir)
 yav=yav/real(nvir)
 zav=zav/real(nvir)
 
 vxav=vxav/real(nvir)
 vyav=vyav/real(nvir)
 vzav=vzav/real(nvir)
 
 vx2av=vx2av/real(nvir)-vxav*vxav
 vy2av=vy2av/real(nvir)-vyav*vyav
 vz2av=vz2av/real(nvir)-vzav*vzav
 
 vrav=vrav/real(nvir)
 vr2av=vr2av/real(nvir)-vrav*vrav

 write(*,*) 'Generated halo with ',ncount,' particles...'
 write(*,*) 'Number of particles inside virial radius: ',nvir
 write(*,*) 'Centre of mass [kpc]: ',xav,yav,zav
 write(*,*) 'Centre of mass velocity [km/s]: ',vxav,vyav,vzav
 write(*,*) 'Velocity dispersion [km/s]: ',sqrt(vx2av),sqrt(vy2av),sqrt(vz2av),sqrt(vx2av+vy2av+vz2av)
 write(*,*) 'Radial velocity (dispersion) [km/s]: ',vrav,sqrt(vr2av)
 
 do i=1,ntype
    nall(i)=0
    highword(i)=0
 end do
 
 if(isbh.eqv..true.) then
    ncount=ncount+1
    np(6)=1
    massarr(6)=mbh
    pos(1,ncount)=0.0
    pos(2,ncount)=0.0
    pos(3,ncount)=0.0
    
    vel(1,ncount)=0.0
    vel(2,ncount)=0.0
    vel(3,ncount)=0.0
 end if
 
 do i=1,ntype
    nall(i)=np(i)
 end do
 
 

#ifdef HDF5    
  if(mode.eq.1.or.mode.eq.2) then
#endif
     call write_gadget_bin(outfile,ncount,nall(1),mode)
#ifdef HDF5  
  else
     call write_hdf5(outfile,ncount,nall(1),mode)
  end if
#endif  
 
end program mk_halo

real(kind=8) function fabel(y)

  use constants
  use params
  use prof_params

  implicit none

  real(kind=8) :: y,p,dp

  integer(kind=4) :: n,nmin,nmax

  nmin=1
  nmax=ntab

  if((ptab(1)-y)*(ptab(nmax)-y).gt.0.0) stop '(ptab(1)-y)*(ptab(nmax)-y).gt.0.0'

  do
     n=(nmin+nmax)/2
     if((ptab(n)-y)*(ptab(nmax)-y).ge.0.0) then
        nmax=n
     else
        nmin=n
     end if

     if(abs(nmax-nmin).le.1) exit
  end do

  dp = drhopsitab(icomp,nmin)*(ptab(nmax)-y)/(ptab(nmax)-ptab(nmin)) + &
       & drhopsitab(icomp,nmax)*(y-ptab(nmin))/(ptab(nmax)-ptab(nmin))

  fabel = dp/dsqrt(xx-y)
    
  return

end function fabel

real(kind=8) function dmass(yy)

    use constants
    use prof_params

    implicit none

    real(kind=8) :: yy

    real(kind=8) :: density,dlr

    integer(kind=4) :: n,nmin,nmax

    nmin=1
    nmax=ntab
    
    if((dlog(rtab(1))-yy)*(dlog(rtab(nmax))-yy).gt.0.0) then
       stop '(rtab(1)-y)*(rtab(nmax)-y).gt.0.0'
    end if
    
    do
       n=(nmin+nmax)/2
       if((dlog(rtab(n))-yy)*(dlog(rtab(nmax))-yy).ge.0.0) then
          nmax=n
       else
          nmin=n
       end if
       
       if(abs(nmax-nmin).le.1) exit
    end do

    dlr=dlog(rtab(nmax))-dlog(rtab(nmin))
    
    density = rhotab(icomp,nmin)*(dlog(rtab(nmax))-yy)/dlr + &
         &        rhotab(icomp,nmax)*(yy-dlog(rtab(nmin)))/dlr

    dmass = 4.*pi*density*dexp(3.*yy)

    return

end function dmass

real(kind=8) function fpot(yy)
  
  use params
  use prof_params
  
  implicit none
  
  real(kind=8) :: y,yy
  real(kind=8) :: dlr,dlru,dlrl,dlpmin,dlpmax
  
  integer(kind=4) :: n,nmin,nmax
  
  y=dexp(yy)

  if(y.gt.rtab(ntab)) then
     fpot = mtab(ntab)/y
     return
  end if
  
  nmin=1
  nmax=ntab
  
  if((dlog(rtab(nmin))-yy)*(dlog(rtab(nmax))-yy).gt.0.0d0) stop '(ptab(1)-y)*(ptab(nmax)-y).gt.0.0'
  
  do
     n=(nmin+nmax)/2
     if((dlog(rtab(n))-yy)*dlog((rtab(nmax))-yy).ge.0.0) then
        nmax=n
     else
        nmin=n
     end if
     
     if(abs(nmax-nmin).le.1) exit
  end do

  dlr=dlog(rtab(nmax))-dlog(rtab(nmin))
  dlru=dlog(rtab(nmax))-yy
  dlrl=yy-dlog(rtab(nmin))
  dlpmin=dlog(mtab(nmin))-2*dlog(rtab(nmin))
  dlpmax=dlog(mtab(nmax))-2*dlog(rtab(nmax))  
  
  fpot = dlpmin*dlru/dlr+dlpmax*dlrl/dlr

  fpot = dexp(fpot)*y
  
  return
  
end function fpot


