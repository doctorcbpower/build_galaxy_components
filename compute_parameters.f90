module compute_parameters
  use structure
contains

  ! Calculate kappa, the epicyclic frequency, given as
  !
  !      kappa^2 = R x dOmega^2/dR + 4 Omega^2
  !
  ! where Omega=vrot/R is the angular velocity. Note that
  !
  ! d/dR (vrot/R)^2 = 1/R^2 dvrot^2/dR + vrot^2 x -2/R^3
  !                 = 2 x vrot/R^2 dvrot/dR -2 x vrot^2/R^3
  !
  ! and so
  !
  !      kappa^2 = 2 x Omega x dvrot/dR -2 x Omega^2 + 4 x Omega^2
  !              = 2 x Omega x (dvrot/dR + Omega)
  
  
  real(kind=8) function get_kappa(r,vrot,dvrotdr)
    implicit none
    real(kind=8) :: r           ! radius
    real(kind=8) :: vrot        ! rotational velocity at that radius
    real(kind=8) :: dvrotdr     ! derivative of rotational velocity at that radius

    get_kappa = 2. * (vrot/r) * (dvrotdr + vrot/r)

    get_kappa = dsqrt(get_kappa)
    
    return
    
  end function get_kappa

  ! Compute the Q parameter for a gas disc; this is defined as
  !
  !     Q = cs*kappa/(pi*G*Sigma)
  !
  
  real(kind=8) function get_qparameter(r,vrot,dvrotdr,cs,sigma)
    implicit none
    real(kind=8) :: r           ! radius
    real(kind=8) :: vrot        ! rotational velocity at that radius
    real(kind=8) :: dvrotdr     ! derivative of rotational velocity at that radius
    real(kind=8) :: cs          ! sound speed
    real(kind=8) :: sigma       ! surface density
    
    get_qparameter=cs*get_kappa(r,vrot,dvrotdr)/pi/ggdt/sigma

    return
    
  end function get_qparameter
  
  real(kind=8) function get_sigma(r)
    implicit none
    real(kind=4) :: r     ! radius
    return
  end function get_sigma

  subroutine get_derivative(x,y,dydx,n)
    implicit none
    
    integer(kind=4) :: n
    real(kind=8), dimension(n) :: x,y,dydx
    
    integer(kind=4) :: i,j
    
    dydx(1)=(y(2)-y(1))/(x(2)-x(1))
    
    do i=2,n-1
       dydx(i)=0.5*((y(i+1)-y(i))/(x(i+1)-x(i))+&
            & (y(i)-y(i-1))/(x(i)-x(i-1)))
    end do
    
    dydx(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
    
    return
    
  end subroutine get_derivative
  
end module compute_parameters
