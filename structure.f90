module structure
  use constants

  !  module physical_properties
    real(kind=8) :: cs      ! sound speed
    real(kind=8) :: rho0    ! central density
    real(kind=8) :: mbh     ! central black hole mass
    real(kind=8) :: radius  ! radius
    real(kind=8) :: m200=0.0  ! virial mass
    real(kind=8) :: mhalo   ! Halo mass
    real(kind=8) :: mvir=0.0  ! virial mass    
    real(kind=8) :: v200=0.0  ! virial velocity
    real(kind=8) :: rs      ! NFW scale radius
    real(kind=8) :: r200    ! virial radius
    real(kind=8) :: c200    ! NFW concentration
    real(kind=8) :: sigma200=0 ! virial velocity dispersion
    real(kind=8) :: sigma0=0.0  ! central surface density
    real(kind=8) :: rdisc   ! disc scale length
    real(kind=8) :: rdisc_ref  ! disc scale reference
    real(kind=8) :: rdisc_hole ! disc scale reference    
    real(kind=8) :: bdisc   ! disc shape parameter (Miyamoto-Nagai disc)
    real(kind=8) :: mbulge  ! bulge mass
    real(kind=8) :: fbulge  ! bulge extent
    real(kind=8) :: rbulge  ! bulge extent  
    real(kind=8) :: mstar   ! stellar mass
    real(kind=8) :: mdisc   ! disc mass
    real(kind=8) :: fdisc=1.0   ! stellar mass
    real(kind=8) :: fhole=0.0   ! size of central disc hole, in units of rdisc
    real(kind=8) :: fdref=0.0   ! size of reference radius, in units of rdisc
    real(kind=8) :: qflat=1.0   ! flattening parameter
    real(kind=8) :: rcore       ! core radius
    real(kind=8) :: rhocrit ! Critical density at z
    real(kind=8) :: temp    ! gas temperature
    real(kind=8) :: fbaryon=0.16 ! baryon fraction
    real(kind=8) :: lambda_b=0.035 ! Bullock spin parameter,
    ! |J|/sqrt(2)/G/M/R
    real(kind=8) :: v_radial=0.0 ! Radial infall velocity

    integer(kind=4) :: prof_type  ! 1=NFW, 2=Hernquist, 3=Burkert, 4=Generalised
    real(kind=8) :: alpha, beta, gamma_slope  ! Exponents of generalised profile
    real(kind=8) :: epsilon,fdecay=3000
    logical :: istaper            ! Taper the mass distribution?

    integer(kind=4), parameter :: ncomp_max=5 ! Maximum number of
    ! N-body components
    integer(kind=4) :: ncomp=1  ! Number of N-body components
    integer(kind=4), dimension(ncomp_max) :: npart_comp=(/0,0,0,0,0/) ! Number of
    ! particles per component
    real(kind=8), dimension(ncomp_max) :: fmass_comp=(/0.,0.,0.,0.,0./) ! Mass fraction
    ! per component, in units of m200
    real(kind=8), dimension(ncomp_max) :: rscale_comp=(/0.,0.,0.,0.,0./) ! Fraction
    ! scale radius per component, in units of rs
    real(kind=8), dimension(ncomp_max) :: raniso_comp=(/0.,0.,0.,0.,0./) ! Fraction
    ! anisotropic radius per component, in units of rs    
    real(kind=8), dimension(ncomp_max) :: rtrunc_comp=(/0.,0.,0.,0.,0./) ! Truncation
    ! radius per component, in units of rvir
    character(kind=1,len=20), dimension(ncomp_max) :: form_comp !
    !Form of component profile; currently NFW or Hernquist

    real(kind=8), dimension(ncomp_max) :: rhostar_comp
    real(kind=8), dimension(ncomp_max) :: beta_comp
    real(kind=8), dimension(ncomp_max) :: rstar_comp

    integer(kind=4), parameter :: nvtab=100
    real(kind=8), dimension(nvtab) :: xtab,ytab,dytab
    real(kind=8), dimension(nvtab) :: vctab,rvctab
    
    !  end module physical_properties

contains
  real(kind=8) function menc(r)
    
    implicit none
    
    real(kind=8) :: r
    real(kind=8) :: fofc

    real(kind=8) :: x,a
    
    menc=0.0
    
#ifdef NFW_HALO
    x=r/rs
    fofc=dlog(1+c200)-c200/(1+c200)
    menc = menc+(mhalo/fofc)*(log(1+x)-x/(1+x))
#endif
    
#ifdef HERNQUIST_HALO
    a=rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
    x=r/a
#ifndef GalIC
    menc=menc+mhalo*(1+a/r200)*(1+a/r200)*x*x/(1+x)/(1+x)
#else
    menc=menc+mhalo*x*x/(1+x)/(1+x)
#endif
#endif
    
#ifdef HERNQUIST_BULGE
    a=fbulge*rs
    x=r/a
    menc=menc+mbulge*x*x/(1+x)/(1+x)      
#endif
    
#ifdef STELLAR_DISC_EXPONENTIAL
    a=rdisc
    x=r/a
    menc=menc+mstar*(1-exp(-x)*(1+x))
#endif      
  end function menc

  real(kind=8) function mhenc(r)
    use nrutils_modules
    implicit none
    
    real(kind=8) :: r
    real(kind=8) :: fofc,fac

    real(kind=8) :: x,a

    if(prof_type.eq.1) then                                ! NFW
       x=r/rs
       fofc=dlog(1+c200)-c200/(1+c200)
       mhenc = (mhalo/fofc)*(dlog(1+x)-x/(1+x))
    else if(prof_type.eq.2) then                          ! Hernquist
       a=rs
#ifdef GalIC       
       a=a*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
#endif
       x=r/a      
       mhenc=mhalo*x*x/(1+x)/(1+x)       
#ifndef GalIC
       mhenc=mhenc*(1+a/r200)*(1+a/r200)
#endif
    else if(prof_type.eq.3) then                          ! Burkert
       x=r/rs
       fofc=dlog(1+c200*c200)+2.*dlog(1+c200)-2.*atan(c200)
       mhenc = (mhalo/fofc)*(dlog(1.+x*x)+2.*dlog(1.+x)-2.*atan(x))
    else
       x=r/rs
       call qromb(integral_gen,dlog(1d-7),dlog(x),fac)       
       mhenc = mhalo * fac
    end if

    return
    
  end function mhenc

  real(kind=8) function mbenc(r)
    implicit none
    
    real(kind=8) :: r,z  ! Project radius, height above midplane
    real(kind=8) :: fofc,phi

    real(kind=8) :: x,a,mtotal
    
#ifdef HERNQUIST_BULGE
    a=fbulge*rs
    x=r/a
    mbenc=mbulge*x*x/(1+x)/(1+x)      
#endif

  end function mbenc
    
  real(kind=8) function mdenc(logr)
    implicit none
    
    real(kind=8) :: logr,r
    
    r=10**logr
    
#ifndef EXPONENTIAL_DISC      
    mdenc=2.*PI*sigma0*(rdisc/r)*r*r/log10e
#else
    mdenc=2.*PI*sigma0*(r*r/(r*r+rdisc_hole*rdisc_hole))*dexp(-(r-rdisc_ref)/rdisc)*r*r/log10e
#endif
  end function mdenc

  real(kind=8) function getrvir(mvir,dvir,rhoc)
    use constants
    implicit none
    real(kind=8) :: mvir,dvir,rhoc
    getrvir=(3*mvir/(4*pi) * (1/dvir) * (1/rhoc))**(1./3.)
    return
  end function getrvir

  real(kind=8) function getmvir(vvir,dvir,rhoc)
    use constants
    implicit none
    real(kind=8) :: vvir,dvir,rhoc
    getmvir=sqrt((3./4./pi) * (1/dvir) * (1/rhoc))
    getmvir=getmvir*(vvir/sqrt(ggdt))**3.0
    return
  end function getmvir

  real(kind=8) function rhodm(r)
    
    implicit none
    
    real(kind=8) :: a,x,r,rho_scale

    if(prof_type.eq.1) then                                ! NFW
        x=r/rs
        rho_scale=mhalo/(dlog(1+c200)-c200/(1+c200))
        rho_scale = rho_scale/(4.*pi*rs**3)
        rhodm = rho_scale/(x*(1+x)*(1+x))
    else if(prof_type.eq.2) then                           ! Hernquist
        a=rs!*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
        x=r/a
#ifndef GalIC
        rho_scale = (mhalo/2./pi/a**3)*((1+r200/a)**2)/(r200/a)**2
#else
        rho_scale = (mhalo/2./pi/a**3)
#endif
        rhodm = rho_scale/(x*(1+x)*(1+x)*(1+x))
     else if(prof_type.eq.3) then                          ! Burkert
        x=r/rs
        rho_scale=mhalo/(dlog(1+c200*c200)+2.*dlog(1+c200)-2.&
             &*atan(c200))
        rho_scale = rho_scale/(pi*rs**3)
        rhodm = rho_scale/((1+x)*(1+x*x))
     else
        rhodm = rho_general(r)
     endif

    return
    
  end function rhodm

  real(kind=8) function rho_general(r)

    use nrutils_modules
    
    implicit none

    real(kind=8) :: a,x,r,rho_scale,fac

    x=r/rs

    call qromb(integral_gen,dlog(1.d-7),dlog(c200),fac)
    
    rho_scale=mhalo/fac

    rho_scale=rho_scale/(4.*pi*rs**3)

    if(r.le.r200) then
       rho_general=rho_scale*x**(-gamma_slope)
       rho_general=rho_general*(1+x**alpha)**((gamma_slope-beta)/alpha)
    else
       rho_general=rho_scale*c200**(-gamma_slope)
       rho_general=rho_general*(1+c200**alpha)**((gamma_slope-beta)/alpha)
       rho_general=rho_general*(x/c200)**epsilon
       rho_general=rho_general*dexp(-(x/c200-1.)/fdecay)
    end if
    
    return

  end function rho_general

  real(kind=8) function integral_gen(x)
    
    implicit none

    real(kind=8) :: x,logint

    if(dexp(x).le.c200) then    
       integral_gen=dexp(x)**(3.-gamma_slope)
       integral_gen=integral_gen*(1+dexp(x)**alpha)**((gamma_slope-beta)/alpha)
       !logint=(3.-gamma_slope)*x+((gamma_slope-beta)/alpha)*dlog(1.+dexp(x)**alpha)
    else
       !logint=-(gamma_slope)*dlog(c200)+((gamma_slope-beta)/alpha)*dlog(1.+c200**alpha)
       !logint=logint+(2.+epsilon)*x-((x/c200)-1.0)/fdecay
       integral_gen=(c200**(-gamma_slope))/(1+c200**alpha)**((beta-gamma_slope)/alpha)
       integral_gen=integral_gen*(rs*dexp(x)/r200)**epsilon
       integral_gen=integral_gen*dexp(-(rs*dexp(x)-r200)/(fdecay*r200))
       integral_gen=integral_gen*dexp(x)**2
    endif

!    integral_gen=dexp(logint)

    return

  end function integral_gen

  real(kind=8) function rho_trunc(r,rscale,rho_scale,rt,rstar&
       &,rho_star)
    
    implicit none
    
    real(kind=8) :: x,y,r,rscale,rho_scale,rt,rstar,rho_star&
         &,lrho_trunc

    if(r.lt.rt) then    
       x=r/rscale
       if(prof_type.eq.1) then                                ! NFW
          rho_trunc = rho_scale/(x*(1+x)*(1+x))
       else if(prof_type.eq.2) then                           !
          !Hernquist
          rho_trunc = rho_scale/(x*(1+x)*(1+x)*(1+x))
       else if(prof_type.eq.3) then                           !
          !Burkert
          rho_trunc = rho_scale/((1+x)*(1+x*x))
       end if
       return
    end if

    y=r/rstar

    if(y.gt.50) then
       rho_trunc=0.0
       return
    end if
    
    x=r/rt

    rho_trunc = rho_star/x/x

    rho_trunc = rho_trunc * dexp(-y)
    
    return
  end function rho_trunc
  
  real(kind=8) function slopedm(r)
    
    implicit none
    
    real(kind=8) :: a,x,r

    if(prof_type.eq.1) then                       ! NFW
       x=r/rs
       slopedm = -(1+3*x)/(1+x)
    else if(prof_type.eq.2) then                  ! Hernquist
       a=rs!*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
       x=r/a
       slopedm = -(1+4*x)/(1+x)
    else                                          ! Burkert
       x=r/rs
       slopedm = -(x/(1+x)+2*x*x/(1+x*x))
    endif
    
    return
    
  end function slopedm

  subroutine taper_mass_profile(mtot,rmax)

    implicit none

    real(kind=8) :: mtot
    real(kind=8) :: mut_tot  ! Untapered total mass
    real(kind=8) :: mut_trunc ! Untapered mass within rtrunc
    real(kind=8) :: rho0
    real(kind=8) :: rtrunc, ftrunc, rmax
    real(kind=8) :: mu
    integer(kind=4) :: i

    prof_type=1
    
    if(form_comp(1).eq.'hernquist') prof_type=2
    if(form_comp(1).eq.'burkert') prof_type=3        
    
    rtrunc=rtrunc_comp(1)*r200
    
    mhalo=m200
    
    if(prof_type.eq.2.and.mhenc(rtrunc).lt.m200) stop 'Check -DGALIC'
    
    mtot=0.0
    
    do i=1,ncomp
       rtrunc=rtrunc_comp(i)*r200
       
       mhalo=fmass_comp(i)*m200

       prof_type=1
       
       if(form_comp(i).eq.'hernquist') prof_type=2
       if(form_comp(i).eq.'burkert') prof_type=3       

       rtrunc=rtrunc_comp(i)*r200

       rstar_comp(i)=-rtrunc/(2.+slopedm(rtrunc))

       rhostar_comp(i)=rhodm(rtrunc)*exp(-2.-slopedm(rtrunc))

       mtot=mtot+mhenc(rtrunc)+4.*pi*rhostar_comp(i)*rtrunc*rtrunc*&
            & rstar_comp(i)*exp(-rtrunc/rstar_comp(i)-exp(-rmax&
            &/rstar_comp(i)))
    end do

  end subroutine taper_mass_profile

  real(kind=8) function get_rho0(mref,rref,rscale)

    use nrutils_modules
    
    implicit none

    real(kind=8) :: x,mref,rref,rscale,fac

    x=rref/rscale        
    if(prof_type.eq.1) then         ! NFW
       get_rho0 = mref/(dlog(1.+x)-x/(1.+x))/4./pi/rscale**3
    else if(prof_type.eq.2) then    ! Hernquist
       get_rho0 = mref*(1.+x)*(1.+x)/x/x/2./pi/rscale**3
    else if(prof_type.eq.3) then    ! Burkert
       get_rho0 = mref/((dlog(1.+x*x)+2.*dlog(1.+x)-2.*atan(x)))/pi&
            &/rscale**3
    else
       call qromb(integral_gen,dlog(1.d-7),dlog(x),fac)
       get_rho0 = mref/fac/4./pi/rscale**3
       
    end if

    return

  end function get_rho0

  real(kind=8) function get_j(r)
    
    implicit none

    real(kind=8) :: r
    
    get_j = lambda_b * sqrt(2.*ggdt*r*mhenc(r))
    
    return

  end function get_j

  real(kind=4) function get_vel_comps(pos,idim)
    
    implicit none

    real(kind=4), dimension(3) :: pos
    integer(kind=4) :: i,idim
    real(kind=4), dimension(3) :: vel    
    real(kind=8) :: r
    real(kind=8) :: vcirc

    r=0.0
    do i=1,2
       r=r+pos(i)**2
    end do

    r=sqrt(r)
    
    vcirc = get_j(r)/r
  
    vel(1)=-vcirc*sin(atan2(pos(2),pos(1)))
    vel(2)=vcirc*cos(atan2(pos(2),pos(1)))
    vel(3)=0.0
    
    get_vel_comps = vel(idim)
    
    return

  end function get_vel_comps

  real(kind=4) function circular_velocity_logarithmic_halo(vscale,rscale,radius)
    
    implicit none

    real(kind=8) :: rscale   ! Scaling radius 
    real(kind=8) :: vscale   ! Scaling velocity 
    real(kind=8) :: radius   ! Radius

    circular_velocity_logarithmic_halo=vscale*radius/sqrt(rscale*rscale+radius*radius)

    return

  end function circular_velocity_logarithmic_halo

  real(kind=4) function get_mean_molecular_weight(X,Y,Z,temperature)

    implicit none

    real(kind=4) :: X,Y,Z   ! Fractions
    real(kind=8) :: temperature  ! Gas temperature
    real(kind=4) :: mu_ion,mu_e

    if(dlog10(temperature).le.4) then
       mu_ion=1./(X/1.+Y/4.)
       mu_e=0.0
       get_mean_molecular_weight=mu_ion
    else 
       mu_ion=1./(X/1.+Y/4.)
       mu_e=1./(X/1.+2.*Y/4.)
       get_mean_molecular_weight=1./(1./mu_ion+1./mu_e)
    end if

    return

  end function get_mean_molecular_weight

  subroutine compute_disc_circular_velocity

    use nrutils_modules 

    implicit none

    integer(kind=4) :: n
    real(kind=8), parameter :: lrmin=-3.,lrmax=3
    real(kind=8) :: logr,dlogr
    real(kind=8) :: y,sum

    dlogr=(lrmax-lrmin)/real(nvtab-1)
    logr=lrmin
    n=0

    ! First we need to tabulate the inner integral
    do
       if(logr-lrmax.gt.1.d-6) exit
       radius=10**logr
       call qromo(inner_disc_integral,lrmin,lrmax,sum,midpnt)
       n=n+1
       xtab(n)=logr
       ytab(n)=sum!radius*sigma0*bessk1(radius/rdisc)!sum
       write(14,*) (10**xtab(n))/rdisc,ytab(n),&
            & radius*sigma0*bessk1(radius/rdisc)
       logr=logr+dlogr
    end do
    
    ! Now we need to compute the outer integral

     dlogr=(1.7+2.5)/real(nvtab-1)    
     logr=-2.
     !    n=0

    do
       if(logr.gt.1.7) exit
       !       n=n+1
       radius=10**logr
       call mqromb(outer_disc_integral,-3.0d0,3.0d0,xtab,ytab,nvtab,sum,1)
!       rvctab(n)=10**logr
!       vctab(n)= 4*pi*ggdt*radius*sum
!       y=radius/2/rdisc
       print *, logr,2*pi*ggdt*radius*sum
!       print *, rvctab(n),vctab(n),(-pi*ggdt*sigma0*radius*(bessi0(y)*bessk1(y)&
!            & - bessi1(y)*bessk0(y)))
!       write(13,'(3e18.8)') rvctab(n),vctab(n),&
!            & (-pi*ggdt*sigma0*radius*(bessi0(y)*bessk1(y)&
!            & - bessi1(y)*bessk0(y)))
       logr=logr+dlogr
    end do
    
    return
    
  end subroutine compute_disc_circular_velocity
  
  real(kind=8) function inner_disc_integral(x)

    use nrutils_modules
    
    implicit none

    real(kind=8) :: r,x,y
    real(kind=8) :: numerator, denominator
    
    r=10**x

#ifdef EXPONENTIAL_DISC    
!    numerator=sigma0*(r*r/(r*r+rdisc_hole*rdisc_hole))*exp(-(r-rdisc_ref)*y)*r
    numerator=sigma0*exp(-r/rdisc)
!    numerator=1.0
#endif
!    denominator=sqrt(r*r-radius**2)
    
    inner_disc_integral=(r*numerator*bessj0(r*radius))*r/0.4342944819

    return

  end function inner_disc_integral

  real(kind=8) function outer_disc_integral(x,xp,yp,np)

    use nrutils_modules     

    implicit none
    
    integer(kind=4) :: n,np,nhi,nlo
    real(kind=8), dimension(np) :: xp,yp
    real(kind=8) :: x,y,dx,x1,xn,h
    real(kind=8) :: u,v,w
!#ifdef CRAP

    n=1

    do
       if(xp(n).ge.x) exit
       n=n+1
    end do

    if(xp(n).eq.x) then
       u=10**xp(n)
       w=yp(n)
       outer_disc_integral = w*u*bessj1(u*radius)*(u/log10e)
       return
    end if
    
    dx=xp(n)-xp(n-1)
    y = yp(n-1)*(xp(n)-x)+yp(n)*(x-xp(n-1))
    y = y/dx

    u=10**x
    w=y

    outer_disc_integral = w*u*bessj1(u*radius)*(u/log10e)
!#endif
!    u=10**x
!    v=u/rdisc
!    outer_disc_integral = (u*bessk1(v)/sqrt(radius*radius-u*u))*(u/log10(exp(1.d0)))
!    print *, n, radius,x,xp(n),yp(n),y,outer_disc_integral
    
    return
    
  end function outer_disc_integral
  
end module structure
