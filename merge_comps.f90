program merge_comp
  use gadget_header
  use io
  use nrutils_modules
  
  implicit none

  ! Process command line arguments
  character(len=132), dimension(20) :: instring ! command line arguments
  character(len=132), dimension(10) :: filename
  character(len=132) :: outfile  
  integer(kind=4) :: i,j,k,l,m,n,nfile
  integer(kind=4), dimension(5) :: iin=[0,0,0,0,0]
  integer(kind=4) :: iout=0
  logical :: fexist
  integer(kind=4) :: ngas,npart
  integer(kind=4), dimension(10) :: npfile=0,npfile_g=0

  real(kind=4), allocatable :: dummy(:,:)
  integer(kind=4), dimension(ntype) :: npall
  integer(kind=4) :: in_mode=1,out_mode=1

  real(kind=4) :: xmin,xmax,ymin,ymax,zmin,zmax
  
  real(kind=8), dimension(ntype) :: massarr_comb=[0,0,0,0,0,0]

  real(kind=4), allocatable :: r2(:)
  integer(kind=4), allocatable :: indx(:)
  real(kind=4) :: r_taper
  integer(kind=4) :: itaper=0
  real(kind=8) :: mtot,mtaper

  real(kind=4) :: f_reduce
  logical :: taper_profile=.false.
  logical :: reduce_stellar_mass=.false.
  logical :: excise_type=.false.  

  integer(kind=4) :: ireduce=0, istar=0, iexcise=0, type_star=-1, type_excise=-1

  logical :: offset=.false.
  integer(kind=4) :: ioffset=0
  character(kind=1,len=132) :: offset_file
  real(kind=4), dimension(3,100) :: dpos,dvel,dangle  
  namelist /position_offset/ dpos
  namelist /angular_offset/ dangle
  namelist /velocity_offset/ dvel
  
  real(kind=4), dimension(3,3) :: rotmat
  real(kind=4), dimension(3) :: pos0,vel0
  logical :: isrotation

  isrotation = .false.
  
  if(command_argument_count().eq.0) then
#ifdef HDF5
     write(*,*) 'Usage: merge_comps.exe -in <infile1> ... -out <outfile>'
     write(*,*) 'I/O options: (-in_hdf5|-out_hdf5)'
     write(*,*) 'I/O options: (-in_swift|-out_swift)'
     write(*,*) 'I/O options: (-out_snap1|-out_snap2)'
#else
     write(*,*) 'Usage: merge_comps.exe -in <infile1> ... -out <outfile>'
     write(*,*) 'I/O options: (-out_snap1|-out_snap2)'
#endif     
     write(*,*) 'Options: (-reduce_stellar_mass <type> -reduce_stellar_mass_fraction <f_s>)'
     write(*,*) 'Options: (-taper_profile -taper_mass <f_m>)'
     write(*,*) 'Options: (-apply_offsets <offset file>)'
     write(*,*) 'Options: (-excise_type <type>)'
     write(*,*) 'Options: (-write_mass_block)'     
     stop
  end if

  n=1
  nfile=0
  
  do
     call get_command_argument(n, instring(n))
     if (len_trim(instring(n)).eq.0) exit  

     select case (instring(n))
       case ("-in")
          nfile=nfile+1
          iin(nfile)=n+1
       case ("-out")
          iout=n+1
#ifdef HDF5
       case ("-in_hdf5")
          in_mode=3
       case ("-out_hdf5")
          out_mode=3
       case ("-in_swift")
          in_mode=4
       case ("-out_swift")
          out_mode=4
#endif
       case ("-out_snap1")
          out_mode=1
       case ("-out_snap2")
          out_mode=2
       case ("-apply_offsets")
          offset=.true.
          ioffset=n+1
       case ("-taper_profile")
          taper_profile=.true.
       case ("-taper_mass")
          itaper=n+1
       case ("-reduce_stellar_mass")
          reduce_stellar_mass=.true.
          istar=n+1
       case ("-reduce_stellar_mass_fraction")
          ireduce=n+1
       case ("-excise_type")
          excise_type=.true.
          iexcise=n+1          
       case ("-write_mass_block")
          ismassblock=.true.
       end select
     n = n+1
  end do

  if(nfile.eq.0) stop 'No input files defined'
  
  n=0
  m=0
  
  do
     if(n.ge.nfile) exit
     n=n+1
     m=m+1
     filename(n)=trim(instring(iin(m)))
     inquire(file=filename(n),exist=fexist)
     if(fexist.eqv..false.) then
        write(*,*) 'Cannot find ',filename(n)
        stop
     end if
  end do

  if(iout.gt.0) then
     outfile=trim(instring(iout))
  else
     if(out_mode.eq.1) then
        outfile='out.snp1.gdt'
     else if(out_mode.eq.2) then
        outfile='out.snp2.gdt'
     else if(out_mode.eq.3) then
        outfile='out.gdt.hdf5'
     else 
        outfile='out.swft.hdf5'        
     end if
  end if

  if(taper_profile) then
     if(itaper.eq.0) stop 'Need to define taper mass'
     read(instring(itaper),*) mtaper
  end if

  if(reduce_stellar_mass) then
     if(istar.eq.0) stop 'Need to define stellar type'
     read(instring(istar),*) type_star
     if(type_star.gt.ntype-1) stop 'star type does not make sense'
     if(ireduce.eq.0) stop 'Need to define stellar mass fraction'
     read(instring(ireduce),*) f_reduce
     if(f_reduce.eq.0.0) then
        reduce_stellar_mass=.false.
        excise_type=.true.
        type_excise=type_star
     end if
  end if

  if(excise_type) then
     if(type_excise.eq.-1.and.iexcise.eq.0) stop 'Need to define excise type'
     if(type_excise.eq.-1) read(instring(iexcise),*) type_excise
     if(type_excise.gt.ntype-1) stop 'excise type does not make sense'
  end if
  
  if(offset) then
     read(instring(ioffset),'(A)') offset_file
     inquire(file=offset_file,exist=fexist)
     if(fexist.eqv..false.) stop 'Could not find offsets file'
     open(1,file=offset_file,status='old')
     rewind(1)
     read(1,nml=position_offset)
     rewind(1)
     read(1,nml=angular_offset)
     rewind(1)     
     read(1,nml=velocity_offset)
     close(1)

!     rotmat(1,1)=cos(dangle(3,1))*cos(dangle(1,1))-cos(dangle(2,1))*sin(dangle(1,1))*sin(dangle(3,1))
!     rotmat(2,1)=cos(dangle(3,1))*sin(dangle(1,1))+cos(dangle(2,1))*cos(dangle(1,1))*sin(dangle(3,1))
!     rotmat(3,1)=sin(dangle(3,1))*sin(dangle(1,1))

!     rotmat(1,2)=-sin(dangle(3,1))*cos(dangle(1,1))-cos(dangle(2,1))*sin(dangle(1,1))*cos(dangle(3,1))
!     rotmat(2,2)=-sin(dangle(3,1))*sin(dangle(1,1))+cos(dangle(2,1))*cos(dangle(1,1))*cos(dangle(3,1))
!     rotmat(3,2)=cos(dangle(3,1))*sin(dangle(2,1))

!     rotmat(1,3)=sin(dangle(2,1))*sin(dangle(1,1))
!     rotmat(2,3)=-sin(dangle(2,1))*cos(dangle(1,1))
!     rotmat(3,3)=cos(dangle(2,1))     

     rotmat(1,1)=1.!cos(dangle(1,1)*pi/180)
     rotmat(2,1)=0.!-cos(dangle(2,1))*cos(dangle(1,1))
     rotmat(3,1)=0.!sin(dangle(1,1)*pi/180)

     rotmat(1,2)=0.0
     rotmat(2,2)=1.0!sin(dangle(2,1))
     rotmat(3,2)=0.0!cos(dangle(2,1))

     rotmat(1,3)=0.!-sin(dangle(1,1)*pi/180)
     rotmat(2,3)=0.0!-cos(dangle(2,1))*sin(dangle(1,1))
     rotmat(3,3)=1.!cos(dangle(1,1)*pi/180)     
  end if
  
  npart=0
  ngas=0

  do i=1,ntype
     npall(i)=0
  end do

  do i=1,n
     ! Read in GADGET data file...
#ifdef HDF5
     if(in_mode.eq.3.or.in_mode.eq.4) then
        call read_hdf5_header(filename(i),npfile(i),npfile_g(i),in_mode)
     else
#endif        
        call read_gadget_bin_header(filename(i),npfile(i),npfile_g(i),in_mode)
#ifdef HDF5
     end if
#endif
     ngas=ngas+npfile_g(i)
     npart=npart+npfile(i)

     do j=1,ntype
        npall(j)=npall(j)+nall(j)
        if(massarr(j).gt.0.0) massarr_comb(j)=massarr(j)
     end do
  end do

  write(*,*) 'Allocating memory for ',npart,' particles...'

  allocate(pos(3,npart))
  if(allocated(pos).eqv..false.) stop 'Could not allocate memory (pos)'
  allocate(vel(3,npart))
  if(allocated(vel).eqv..false.) stop 'Could not allocate memory (vel)'  
  allocate(id(npart))
  if(allocated(id).eqv..false.) stop 'Could not allocate memory (id)'  
  allocate(ptype(npart))
  if(allocated(ptype).eqv..false.) stop 'Could not allocate memory (ptype)'  
  allocate(mass(npart))
  if(allocated(mass).eqv..false.) stop 'Could not allocate memory (mass)'  
  if(ngas.gt.0) then
     allocate(ug(ngas))
     if(allocated(ug).eqv..false.) stop 'Could not allocate memory (ug)'
     allocate(hsml(ngas))
     if(allocated(hsml).eqv..false.) stop 'Could not allocate memory (hsml)'
  end if
  
  k=0
  l=0

  do i=1,n
#ifdef HDF5
     if(in_mode.le.2) then
#endif
        call read_gadget_bin(filename(i),npfile(i),npfile_g(i),k,in_mode)
#ifdef HDF5
     else
        call read_hdf5(filename(i),npfile(i),npfile_g(i),k,in_mode)
     end if
#endif
     if(offset) then
        if(i.ge.1) then
           do j=1,npfile(i)
              do m=1,3
                 pos0(m)=pos(m,j+l)
                 vel0(m)=vel(m,j+l)
              end do
              
              pos(1,j+l)=rotmat(1,1)*pos0(1)+rotmat(2,1)*pos0(2)+rotmat(3,1)*pos0(3)
              pos(2,j+l)=rotmat(1,2)*pos0(1)+rotmat(2,2)*pos0(2)+rotmat(3,2)*pos0(3)
              pos(3,j+l)=rotmat(1,3)*pos0(1)+rotmat(2,3)*pos0(2)+rotmat(3,3)*pos0(3)

              vel(1,j+l)=rotmat(1,1)*vel0(1)+rotmat(2,1)*vel0(2)+rotmat(3,1)*vel0(3)
              vel(2,j+l)=rotmat(1,2)*vel0(1)+rotmat(2,2)*vel0(2)+rotmat(3,2)*vel0(3)
              vel(3,j+l)=rotmat(1,3)*vel0(1)+rotmat(2,3)*vel0(2)+rotmat(3,3)*vel0(3)                               

              do m=1,3
                 pos(m,j+l)=pos(m,j+l)+dpos(m,i)
                 vel(m,j+l)=vel(m,j+l)+dvel(m,i)
              end do
           end do
        end if
     end if
     l=l+npfile(i)

!     k=k+npfile(i)
  end do
  
  do j=1,ntype
     massarr(j)=massarr_comb(j)
  end do
  
  call sort_particles_by_ptype(npart,ngas,npall)
     
  write(*,*) 'Writing data to file'

  if(out_mode.eq.4) then
     xmin=1.e10
     ymin=1.e10
     zmin=1.e10
     xmax=-1.e10
     ymax=-1.e10
     zmax=-1.e10     
     do i=1,npart
        xmin=min(xmin,pos(1,i))
        xmax=max(xmax,pos(1,i))
        ymin=min(ymin,pos(2,i))
        ymax=max(ymax,pos(2,i))        
        zmin=min(zmin,pos(3,i))
        zmax=max(zmax,pos(3,i))
     end do
     BoxSize=1.1*max(xmax-xmin,ymax-ymin,zmax-zmin)
     do i=1,npart
        do j=1,3
           pos(j,i)=pos(j,i)+0.5*BoxSize
        end do
     end do
  end if

  if(excise_type) then
     l=0
     m=0
     do i=1,ntype
        if(i.ne.type_excise) then
           do j=l+1,l+np(i)
              m=m+1
              do n=1,3
                 pos(m,n)=pos(j,n)
                 vel(m,n)=vel(j,n)
              end do
              mass(m)=mass(j)
              id(m)=id(j)
              ug(m)=ug(j)
           end do
        end if
        l=l+np(i)
        if(i.eq.type_excise) then
           npart=npart-nall(i)
           np(i)=0
           nall(i)=0
           massarr(i)=0.0
        end if
     end do
  end if

  if(taper_profile) then
     allocate(indx(npart),r2(npart))

     mtot=0.0
     do i=1,npart
        mtot=mtot+mass(i)
        r2(i)=pos(1,i)**2+pos(2,i)**2+pos(3,i)**2
     end do
     
     call indexx(npart,r2,indx)
     
     r_taper=0.0
     
     i=1
     mtaper=mtaper*mtot
     mtot=0.0
     
     do
        j=indx(i)
        mtot=mtot+mass(j)
        if(mtot.gt.mtaper) exit
        i=i+1
     end do
     
     j=i-1
     r_taper=sqrt(r2(indx(j)))
     
     k=187619
     n=0
     do i=1,npart
        l=indx(i)
        if(ran3(k).lt.exp(-(sqrt(r2(l))-r_taper)/r_taper)) then
           n=n+1
           indx(n)=l
        end if
     end do

     npart=n
     
     call sample_particles_by_id(npart,npall,indx)     

     deallocate(indx,r2)
  end if

  if(reduce_stellar_mass) then
     
     i=1
     j=0
     k=0
     do
        if(i.ge.type_star+1) exit
        j=k+1
        k=j+np(i)
        i=i+1
     end do

     do i=j,k
        mass(i)=f_reduce*mass(i)
     end do

     massarr(type_star)=f_reduce*massarr(type_star)
  end if

#ifdef HDF5  
  if(out_mode.eq.3.or.out_mode.eq.4) then
     call write_hdf5(outfile,npart,npall(1),out_mode)
  else
#endif     
     call write_gadget_bin(outfile,npart,npall(1),out_mode)  
#ifdef HDF5
  end if
#endif  

  deallocate(pos,vel,id,mass,ptype)
  if(ngas.gt.0) deallocate(ug)

end program merge_comp

