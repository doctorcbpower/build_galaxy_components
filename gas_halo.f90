module gas_halo
  use nrutils_modules
  use constants
  use structure
  
  real(kind=8) :: a0       ! Normalisation constant
  real(kind=8) :: fac
  real(kind=8) :: p_index  ! Polytropic index of halo gas
  real(kind=8) :: conc     ! Concentration
contains
  real(kind=8) function mhenc_norm(x)
    implicit none
    
    real(kind=8) :: x
    
    if(prof_type.eq.1) then
       mhenc_norm = (log(1+x)-x/(1+x))
    else if(prof_type.eq.2) then
       mhenc_norm = x*x/(1+x)/(1+x)
    end if
    
    return
  end function mhenc_norm

  real(kind=8) function rhodm_norm(x)
    implicit none
    
    real(kind=8) :: x
    
    if(prof_type.eq.1) then
       rhodm_norm = 1./x/(1+x)/(1+x)
    else if(prof_type.eq.2) then
       rhodm_norm = 1./x/(1+x)/(1+x)/(1+x)
    end if
    
    return
  end function rhodm_norm
  
  real(kind=8) function mhenc_over_r2(logx)
    implicit none
    
    real(kind=8) :: logx,y

    y=dexp(logx)

    mhenc_over_r2 = dexp(dlog(mhenc_norm(y))-2*dlog(y))*dexp(logx)

    return

  end function mhenc_over_r2
  
  real(kind=8) function func(x)
    
    implicit none
    
    real(kind=8) :: x
    
    func = x*x*rhog(x) 
    
    return
  end function func
  
  real(kind=8) function rhog(x)
    
    implicit none
    
    real(kind=8) :: x
    real(kind=8), parameter :: xmin=1d-6
    real(kind=8) :: sum

    call qromb(mhenc_over_r2,dlog(xmin),dlog(x),sum)

    rhog = (1.0 - (3./a0) * fac * sum)**(1./(p_index-1.))

    return
    
  end function rhog

  ! Determine the mass profile of a gaseous halo in hydrostatic equilibrium, 
  ! embedded in a dark matter halo.
  !
  ! Assume that the gas has a polytropic equation of state, P=A x rho^gamma, and
  ! that the gas has a density profile rho(r) = rho0 f(r), which is different
  ! to the underlying dark matter density profile.
  !
  ! In this case, the condition for hydrostatic equilibrium can be written as
  !
  !        1/rho x dP/dr = -G M(<r)/r^2
  !
  ! and so, assuming an isothermal equation of state, we get
  !
  !        dlog(rho(r))/dr = -G mu m_p/(k T_g) x M(r)/r^2
  !
  ! which gives 
  !
  !        rho(r) = rho0xexp(fac x Int)
  !
  ! where
  !
  !        fac = -G mu m_p/(k T_g) = -3.2775 (mu/0.63) (T_g/10^6K) x
  !                                        (Mvir/10^10 Msol) (rvir/kpc)^-1 x
  !                                                c/g(c)
  !
  ! for the case of a NFW halo.
  !
  ! We need to evaluate the constant A and we do this by requiring that the
  ! gas and dark matter density profiles have matching logarithmic slopes at some
  ! fiducial radius -- the virial radius, for example. The requirement is that 
  ! the radius be large and well outside the central region of the halo.
  ! This allows the constant A to be fixed.
  !
  ! In this case
  !
  !       d log(rho)/alog(r) = r x (1/gamma-1) x 
  !              [1 - GxAx(gamma-1)/gamma int^R_0 M(<r)/r^2 dr]^-1  x
  !               GxAx(gamma-1)/gamma x d/dr (int^R_0 M(<r)/r^2 dr);
  !
  ! we can make use of Leibnitz's rule to evaluate the derivative of the intergal
  ! int^R_0 M(<r)/r^2 dr. This gives M(<r)/r^2, and so
  !
  !     d log(rho)/alog(r) = (1/gamma-1) x 
  !              [1 - GxAx(gamma-1)/gamma int^R_0 M(<r)/r^2 dr]^-1  x
  !               GxAx(gamma-1)/gamma x M(<r)/r.
  ! 
  ! Matching this to the dark matter slope S gives
  !
  !     A =               S x gamma
  !        ------------------------------------------------------
  !        G x M(<r)/r + x G x (gamma-1) int^R_0 M(<r)/r^2 dr x S
  ! 
  ! If the two profiles are to agree over a larger range of radii, this requires
  ! that A be a shallow function of gamma -- effectively flat.
  !
  !
  ! Chris Power, 10th February 2009
  
  subroutine menc_gas(xmin,xmax,mgas,f_g,ntab,rtab,mtab,utab)
    implicit none
    
    real(kind=8) :: xmin,xmax    ! Normalised radial bounds, in units of r200
    real(kind=8) :: x,x1,x2,xfac ! Normalised radial bounds, in units of rs
    real(kind=8) :: f_g          ! Baryon fraction
    real(kind=8) :: mgas         ! Gas mass
    real(kind=8) :: temp0        ! Central gas temperature
    real(kind=8) :: f1,f2        ! Multiplicative factors 
    real(kind=8) :: rho0g        ! Central gas density
    real(kind=8) :: rho0dm       ! Central dark matter density
    real(kind=8) :: rscale       ! Generic scale radius 
    real(kind=8) :: sum,sum1
    
    integer(kind=4) :: ntab
    real(kind=8), dimension(ntab) :: rtab,mtab,utab
    
    integer(kind=8) :: i

    molecular_weight = 0.61
    
    if(xmin.lt.0) stop 'xmin.lt.0'        ! Remember, xmin and xmin are in units
                                          ! of r200
    if(xmax.lt.0) stop 'xmax.lt.0'
    
    if(xmin.ge.xmax) stop 'xmin.ge.xmax'

    rscale=rs

!    if(prof_type.eq.2) rscale=rscale*sqrt(2*(log(1+c200)-c200/(1+c200)))
    
    conc=r200/rscale
    
    ! First we need to evaluate the constant eta0 (eq 23 in Komatsu & Seljak 2001); 
    !     eta0 = 1/gamma x [ (-3/S*) x (f(x*)/x*)/(f(c)/c) + 
    !                            3x(gamma-1)x(c/f(c)xint^x*_0 f(u)/u**2 du]
    ! where
    !
    ! * gamma is the polytropic index - set by requiring DM and gas density profile slopes
    !   match at large radius, which we assume is R200
    ! * S* is the slope of the DM density profile at R200
    ! * x* is this large normalising radius (set to R200) in units of the scale radius
    ! * f(x) is the normalised enclosed mass profile - which is defined in mhenc_norm()
    !

    p_index=1.15 + 0.01*(conc-6.5) ! Eq. 25 from Komatsu & Seljak (2001)
    
    print *, 'Using a gamma of ',p_index
    
    x1=dlog(conc)+dlog(xmin)   ! These have been passed in in units of r200, 
    x2=dlog(1.*conc)           ! so we're now expressing in units of rs

    call qromo(mhenc_over_r2,x1,x2,sum,trapzd)
    
    a0 = 3.*(p_index-1)*(conc/mhenc_norm(conc))*sum
    a0 = a0+(-3./slopedm(r200))*(mhenc_norm(dexp(x2))/dexp(x2))*(conc/mhenc_norm(conc))
    a0 = (1.0/p_index)*a0
    
    ! From equation 20 of Komatsu & Seljak 2001, we can work out the central temperature.
    !
    !    1./eta0 = G * mu * m_p * M200 / (3 * r200 * kb * T0)
    !
    ! and so
    !
    !       T0 = (G * mu * m_p * M200)/(3 * r200 *kb) * eta0
    !
    ! Note that this can be recast as
    !
    !      kT0 = 1/3 * mu*mp * V200**2
    !
    ! Make sure that V200 is expressed in m/s rather than km/s for units to make sense!
    
    temp0 = (1./3.) * (molecular_weight*mp/kb) * (kms_in_ms*kms_in_ms) * (ggdt*m200/r200) * a0

    write(*,*) 'Central halo temperature [K]: ',temp0
    
    ! Now we need to evaluate the central density rho0 by integrating
    !
    ! rho(r) = rho0 x [1 - (3/eta0) x (gamma-1)/gamma x
    !                     c/f(c) x int^r/rs_0 M(<u)/u^2 du ]^(1/gamma-1)
    !
    ! over the interval (0,rvir).
    !
    ! This gives us
    !
    !     f_g x Mvir = 4 x PI x rho0 x rs^3 x int_0^c rho(u) u^2 du
    !

    fac=(1.-1./p_index) * (conc/mhenc_norm(conc))

    call qromb(func,dexp(x1),dexp(x2),sum)

    rho0g = f_g*m200/(4*pi*(rscale**3))/sum
    
    mgas = 4.*pi*(rscale**3.)*rho0g*sum

    ! Now generate table of radii and masses
    
    xfac = (xmax/1.01/xmin)**(1./real(ntab-1))
    
    x = 1.01*xmin*conc 
    
    do i=1,ntab
       call qromb(func,dexp(x1),x,sum)
       mtab(i) = 4*pi*(rscale**3)*rho0g*sum/mgas
       rtab(i) = x*rscale
       utab(i) = temp0*rhog(x)**(p_index-1)
       utab(i) = utab(i)/u_to_temp/molecular_weight
#ifdef DEBUG
       write(32,*) rtab(i),rtab(i)/r200,rho0g*rhog(x),mtab(i),utab(i)
#endif       
       x = x*xfac
    end do

    return
    
  end subroutine menc_gas

end module gas_halo
