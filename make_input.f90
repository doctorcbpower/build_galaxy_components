module make_input
  use nrutils_modules
contains

  subroutine get_nglass(glassfile,nglass)
    implicit none

    character(kind=1,len=*) :: glassfile
    integer(kind=4), dimension(6) :: np
    integer(kind=4) :: i,nglass

    open(1,file=glassfile,status='old',form='unformatted')  ! It's in SnapFormat=1    
    read(1) np
    close(1)

    do i=1,6
       if(np(i).gt.0) nglass=np(i)
    end do

    return
  end subroutine get_nglass

  subroutine get_glass_particles(glassfile,nvert,nbox,nglass,ntot,xglass,yglass,zglass,mode)
    implicit none

    character(kind=1,len=*) :: glassfile
    integer(kind=4) :: mode
    integer(kind=4) :: ntot,nglass,nside,nvert,nbox

    real(kind=4), dimension(ntot) :: xglass,yglass,zglass
    integer(kind=4), allocatable :: indx(:)
    real(kind=4), allocatable :: r2(:),pos(:,:)

    real(kind=8) :: xcm,ycm,zcm  ! Centre of mass
    real(kind=4) :: lbox,xmax,xmin,ymax,ymin,zmin,zmax  ! Boundaries in (x,y,z) dimension

    integer(kind=4) :: i,j,k,l,m,n

    ! First open up the file to read
    open(1,file=glassfile,status='old',form='unformatted')  ! It's in SnapFormat=1
    read(1)
    allocate(pos(3,nglass))
    read(1) ((pos(i,j),i=1,3),j=1,nglass)
    close(1)

    ! Now increase the number of particles in the shell.

    if(mode.eq.0) then            ! If this is a spherical shell...
       if(nbox.gt.1) then
          n=0
          
          do k=0,nbox-1
             do j=0,nbox-1
                do i=0,nbox-1
                   do m=1,nglass
                      n=n+1
                      xglass(n)=pos(1,m)+real(i)
                      yglass(n)=pos(2,m)+real(j)
                      zglass(n)=pos(3,m)+real(k)
                   end do
                end do
             end do
          end do
          
          nglass=n

          do i=1,nglass
             xglass(i)=xglass(i)/real(nbox)
             yglass(i)=yglass(i)/real(nbox)
             zglass(i)=zglass(i)/real(nbox)
          end do
       else
          do i=1,nglass
             xglass(i)=pos(1,i)
             yglass(i)=pos(2,i)
             zglass(i)=pos(3,i)
          end do          
       end if

       xmin=xglass(1)
       xmax=xglass(1)
       
       ymin=yglass(1)
       ymax=yglass(1)
       
       zmin=zglass(1)
       zmax=zglass(1)
       
       do i=1,nglass
          xmin=min(xmin,xglass(i))
          xmax=max(xmax,xglass(i))
          ymin=min(ymin,yglass(i))
          ymax=max(ymax,yglass(i))
          zmin=min(zmin,zglass(i))
          zmax=max(zmax,zglass(i))
       end do
       
       lbox=max(xmax-xmin,ymax-ymin,zmax-zmin)

       xcm=0.5*(xmin+xmax)
       ycm=0.5*(ymin+ymax)
       zcm=0.5*(zmin+zmax)

       deallocate(pos)

       allocate(pos(3,nglass))
       
       do i=1,nglass
          pos(1,i)=(xglass(i)-xcm)/lbox+0.5
          if(pos(1,i).gt.1.0) pos(1,i)=pos(1,i)-1.0
          if(pos(1,i).lt.0.0) pos(1,i)=pos(1,i)+1.0
          
          pos(2,i)=(yglass(i)-ycm)/lbox+0.5
          if(pos(2,i).gt.1.0) pos(2,i)=pos(2,i)-1.0
          if(pos(2,i).lt.0.0) pos(2,i)=pos(2,i)+1.0

          pos(3,i)=(zglass(i)-zcm)/lbox+0.5
          if(pos(3,i).gt.1.0) pos(3,i)=pos(3,i)-1.0
          if(pos(3,i).lt.0.0) pos(3,i)=pos(3,i)+1.0
       end do
       
!!! Now sort particles with respect to the centre of the cube ...
       
       allocate(r2(nglass),indx(nglass))
       
       do i=1,nglass
          r2(i)=(pos(1,i)-0.5)**2.0+(pos(2,i)-0.5)**2.0+(pos(3,i)-0.5)**2.0
       end do
       
       call indexx(nglass,r2,indx)

       i=1
       
       nglass=0
       
       do
          if(r2(indx(i)).gt.(0.5*0.5)) exit
          nglass=nglass+1
          i=i+1
       end do
       
       if(nglass.eq.0) stop 'Warning! No glass particles retained!'
              
       do i=1,nglass
          j=indx(i)
          xglass(i)=pos(1,j)
          yglass(i)=pos(2,j)
          zglass(i)=pos(3,j)
       end do

       ntot=nglass

       deallocate(pos,indx,r2)
       
       return
    end if
    
    nside=floor(nglass**(1./3.))        ! If this is a disc...
    
#ifdef DEBUG
    write(*,*) 'Number of particles on a side:',nside
#endif
  
    ! Cube is periodic in x,y,z; break symmetry by shaving in z direction.
    
    n=0
    
    allocate(indx(nglass))
    
    do i=1,nglass
       if(pos(3,i).lt.real(nvert)/real(nside)) then
          n=n+1
          indx(n)=i
       end if
    end do
    
    nglass=n
    
    xcm=0.0
    ycm=0.0
    zcm=0.0
    
    xmax=-1.e10
    ymax=-1.e10
    zmax=-1.e10
    xmin=1.e10
    ymin=1.e10
    zmin=1.e10
    
    n=0
    
    do i=0,nbox
       do j=0,nbox
          do k=1,nglass
             n=n+1
             m=indx(k)
             xglass(n)=pos(1,m)+real(i)
             yglass(n)=pos(2,m)+real(j)
             zglass(n)=pos(3,m)
             
             xcm=xcm+xglass(n)
             ycm=ycm+yglass(n)
             zcm=zcm+zglass(n)
             
             xmax=max(xmax,xglass(n))
             ymax=max(ymax,yglass(n))
             zmax=max(zmax,zglass(n))
             xmin=min(xmin,xglass(n))
             ymin=min(ymin,yglass(n))
             zmin=min(zmin,zglass(n))
          end do
       end do
    end do
    
    nglass=n
    
    deallocate(indx)
    deallocate(pos) ! Free up memory for later...
    
    lbox=max(xmax-xmin,ymax-ymin)
    if(abs(real(1+nbox)-lbox)/real(1+nbox).gt.1.e-2) stop 
    
    xcm=xcm/real(nglass)
    ycm=ycm/real(nglass)
    zcm=zcm/real(nglass)
    
#ifdef DEBUG
    write(*,*) 'Read ',nglass,' glass particles...'
    write(*,*) 'Centre of mass ',xcm,ycm,zcm
#endif
    
    allocate(r2(nglass),pos(3,nglass))
    
    do i=1,nglass
       pos(1,i)=2.*(xglass(i)-xcm)/real(1+nbox)
       pos(2,i)=2.*(yglass(i)-ycm)/real(1+nbox)
       pos(3,i)=(zglass(i)-zcm)
       r2(i)=pos(1,i)**2+pos(2,i)**2
    end do
    
    allocate(indx(nglass))
    
    call indexx(nglass,r2,indx)
    
    i=0
    
    do
       i=i+1
       if(r2(indx(i)).gt.1) exit
    end do

    nglass=i

    xmin=1.e10;xmax=-1.e10;ymin=1.0e10;ymax=-1.0e10
    do i=1,nglass
       j=indx(i)
       xglass(i)=pos(1,j)/sqrt(r2(j))
       yglass(i)=pos(2,j)/sqrt(r2(j))
       zglass(i)=pos(3,j)/(0.5*real(nvert)/real(nside))
       xmin=min(xglass(i),xmin)
       xmax=max(xglass(i),xmax)       
       ymin=min(yglass(i),ymin)
       ymax=max(yglass(i),ymax)       

       if(abs(zglass(i)).gt.1.0) then
          if(zglass(i).gt.1.0) zglass(i)=zglass(i)-1.0
          if(zglass(i).lt.-1.0) zglass(i)=zglass(i)+1.0          
       endif
    end do

    print *, xmin,xmax
    print *,ymin,ymax

#ifdef DEBUG
    open(33,file='input_glass_particles.txt',status='unknown')
    do i=1,nglass
       write(33,*) xglass(i),yglass(i),zglass(i)
    end do
    close(33)
#endif

    deallocate(r2,indx,pos)
    
    write(*,*) 'Using ',nglass,' glass particles...'

    ntot=nglass
    
    return
  end subroutine get_glass_particles
end module make_input
