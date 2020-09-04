module test
  integer(kind=4) :: nrbin,nzbin
  real(kind=8), allocatable :: rtab(:),ztab(:)
  real(kind=8), allocatable :: phistar(:,:),pot(:,:)
end module test

subroutine compute_background_potential
  use nrutils_modules
  use constants
  use structure
  use test
  implicit none
  
  real(kind=8) :: z  ! Projected radius, height above midplane
  real(kind=8) :: fofc,phi,lphi
  
  real(kind=8) :: x,a,mtotal
  
  integer(kind=4) :: i,j
  
  real(kind=8), parameter :: bmin=0.0,bmax=20.
  
  real(kind=8), external :: dint
  
  real(kind=8) :: sum
  
  allocate(pot(nrbin,nzbin))
  
  do i=1,nrbin
     do j=1,nzbin
        pot(i,j)=0.0
     end do
  end do

#ifdef LOGARITHMIC_HALO
  do i=1,nrbin
     do j=1,nzbin
        phi = 0.5*v200*v200*dlog(rcore*rcore+rtab(i)*rtab(i)+(ztab(j)*ztab(j)/qflat/qflat))
        pot(i,j)=pot(i,j)+phi
     end do
  end do
  write(*,*) 'Finished constructing potential of Logarithmic Halo'      
#endif
  
#ifdef NFW_HALO
  fofc=dlog(1+c200)-c200/(1+c200)
  
  do i=1,nrbin
     do j=1,nzbin
        x=dsqrt(rtab(i)*rtab(i)+ztab(j)*ztab(j))/rs
        lphi = dlog(ggdt*mhalo/(rs*fofc))+dlog(dlog(1.0+x))-dlog(x)
        phi = -dexp(lphi)
        pot(i,j)=pot(i,j)+phi
     end do
  end do
  write(*,*) 'Finished constructing potential of NFW Halo'      
#endif
  
#ifdef HERNQUIST_HALO
  a=rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
#ifndef GalIC
  mtotal=mhalo*(1+a/r200)*(1+a/r200)
#else
  mtotal=mhalo
#endif

  do i=1,nrbin
     do j=1,nzbin
        x=dsqrt(rtab(i)*rtab(i)+ztab(j)*ztab(j))/a
        lphi=dlog(ggdt*mtotal/a)-dlog(1+x)
        phi =-dexp(lphi)
        pot(i,j)=pot(i,j)+phi
     end do
#ifdef WRITE_ESCAPE_VELOCITY         
     write(33,*) rtab(i)/r200,sqrt(2*abs(pot(i,1)))
#endif         
  end do
  write(*,*) 'Finished constructing potential of Hernquist Halo'
#endif
  
#ifdef HERNQUIST_BULGE
  a=rbulge
  mtotal=mbulge
  
  do i=1,nrbin
     do j=1,nzbin
        x=dsqrt(rtab(i)*rtab(i)+ztab(j)*ztab(j))/a
        lphi=dlog(ggdt*mtotal/a)-dlog(1+x)
        phi = -dexp(lphi)
        pot(i,j)=pot(i,j)+phi
     end do
#ifdef WRITE_ESCAPE_VELOCITY         
     write(34,*) rtab(i)/r200,sqrt(2*abs(pot(i,1)))         
#endif
  end do
  write(*,*) 'Finished constructing potential of Hernquist Bulge'      
#endif
  
#ifdef STELLAR_DISC_EXPONENTIAL
  a=rdisc
  mtotal=mstar
  
  do i=1,nrbin
     radius=rtab(i)
     do j=1,nzbin
        call qqromb(dint,bmin,bmax,ztab(j),sum)
        pot(i,j)=pot(i,j)-ggdt*mtotal*sum            
     end do
#ifdef WRITE_ESCAPE_VELOCITY                  
     write(35,*) rtab(i)/r200,sqrt(2*abs(pot(i,1)))
#endif
  end do
  write(*,*) 'Finished constructing potential of Exponential Stellar Disc'      
#endif
  
  return
  
end subroutine compute_background_potential


