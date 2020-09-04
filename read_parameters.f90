module read_params
  use structure
contains
  subroutine read_parameters(infile,ndisc,mbh,outfile)
    implicit none
    
    character(kind=1,len=*) :: infile,outfile
    integer(kind=4) :: ndisc
    real(kind=4) :: mbh

    integer(kind=4) :: n
    integer(kind=4), parameter :: nmax=100
    character(kind=1,len=10) :: instring
    real(kind=4) :: var

    print *, infile
    
    n=0
    
    open(1,file=infile,status='old')

    do
       if(n.gt.nmax) stop 'Error: n.gt.nmax'
       read(1,'(a10,f18.8)',end=10) instring,var
       

!!$    read(1,*) m200
!!$    read(1,*) c200
!!$    read(1,*) mdisc
!!$    read(1,*) fdisc
!!$    read(1,*) mbulge
!!$    read(1,*) fbulge
!!$    read(1,*) mstar
!!$    read(1,*) temp
!!$    read(1,*) ndisc
!!$    read(1,*) mbh
!!$    read(1,'(a)') outfile
    close(1)

10  continue
    
    stop
    
  end subroutine read_parameters
end module read_params
