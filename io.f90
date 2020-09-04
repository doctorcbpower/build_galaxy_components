module io
  use gadget_header
  use constants
#ifdef HDF5  
  use hdf5
#endif
  use nrutils_modules
  real(kind=4), allocatable :: pos(:,:)   ! Positions
  real(kind=4), allocatable :: vel(:,:)   ! Velocities
  integer(kind=4), allocatable :: id(:)   ! Particle IDs
  integer(kind=4), allocatable :: ptype(:) ! Particle Type
  real(kind=4), allocatable :: mass(:)    ! Masses
  real(kind=4), allocatable :: ug(:)      ! Internal Energies
  real(kind=4), allocatable :: hsml(:)    ! Smoothing Lengths
  logical:: ismassblock
contains
  
  subroutine write_gadget_bin(outfile,npart,ngas,mode)
    implicit none

    character(kind=1,len=*) :: outfile     ! output file in GADGET binary format
    integer(kind=4) :: npart,ngas          ! number of particles, gas particles
    integer(kind=4) :: mode         ! mode=1 for SnapFormat=1, mode=2 for SnapFormat=2
    integer(kind=4) :: write_masses        ! write to block if 1
    integer(kind=4) :: i,j,k               ! Dummy variables
    character(kind=1,len=4) :: tag         ! Descriptor for GADGET SnapFormat 2

1234 format(a4)

    j=0
    k=0
    
    do i=1,ntype
       if(np(i).gt.0) then
          k=k+1
          if(massarr(i).gt.0.0) j=j+1
       end if
    end do

    write_masses=1

    if(j.eq.k) then
       write_masses=0
    else
       j=0
       do i=1,ntype
          if(np(i).gt.0.and.massarr(i).gt.0.0) then
             do k=j+1,j+np(i)
                mass(k)=massarr(i)
             end do
             massarr(i)=0.0
          end if
          j=np(i)
       end do       
    end if
    
    write(*,*) 'Writing GADGET binary file to ',trim(outfile)
    
    open(1,file=outfile,status='unknown',form='unformatted')

    write(tag,1234) 'HEAD'
    if(mode.eq.2) write(1) tag,256
    write(1) np,massarr,time,1./time-1.,flagsfr,flagfeedback,nall,&
         &     flagcooling,NumFiles,BoxSize,Omega0,OmegaLambda,HubbleParam,&
         &     flagstellarage,flagmetals,highword,flagentropy,flagdblprc,&
         &     unused

    write(tag,1234) 'POS'    
    if(mode.eq.2) write(1) tag,3*npart*4
    write(*,*) 'Writing positions...'
    write(1) ((pos(j,i),j=1,3),i=1,npart)

    write(tag,1234) 'VEL'    
    if(mode.eq.2) write(1) tag,3*npart*4
    write(*,*) 'Writing velocities...'
    write(1) ((vel(j,i),j=1,3),i=1,npart)

    write(tag,1234) 'IDS'        
    if(mode.eq.2) write(1) tag,npart*4
    write(*,*) 'Writing identities...'
    write(1) (id(i),i=1,npart)

    if(write_masses.eq.1) then
       write(tag,1234) 'MASS'    
       if(mode.eq.2) write(1) tag,npart*4
       write(*,*) 'Writing masses...'
       write(1) (mass(i),i=1,npart)
    end if

    if(ngas.gt.0) then
       write(tag,1234) 'U'           
       if(mode.eq.2) write(1) tag,np(1)*4
       write(*,*) 'Writing internal energies...'
       write(1) (ug(i),i=1,ngas)
    end if

    close(1)
    write(*,*) 'Finished writing output'
    return
  end subroutine write_gadget_bin

  subroutine read_gadget_bin_header(infile,npart,ngas,mode)
    implicit none

    character(kind=1,len=*) :: infile     ! output file in GADGET binary format
    integer(kind=4) :: npart,ngas         ! number of particles, gas particles
    integer(kind=4) :: mode               ! mode=1 for SnapFormat=1, mode=2 for SnapFormat=2
    character(kind=1,len=4) :: check_block
    integer(kind=4) :: i                  ! Dummy variable

    mode=1
    open(1,file=infile,status='old',form='unformatted')
    read(1) check_block
    if(check_block.eq.'HEAD') then
       write(*,*) 'Assuming snapshot format 2'
       mode=2
    else
       rewind 1
    end if
    
    read(1) np,massarr,time,time,flagsfr,flagfeedback,nall,&
         &     flagcooling,NumFiles,BoxSize,Omega0,OmegaLambda,HubbleParam,&
         &     flagstellarage,flagmetals,highword,flagentropy,flagdblprc,&
         &     unused

    redshift=time
    expansion=1./(1.+redshift)
    time=expansion
    
    npart=0
    
    do i=1,ntype
       npart=npart+nall(i)
    end do
    
    ngas=np(1)

    close(1)

    return
    
  end subroutine read_gadget_bin_header
    
  subroutine read_gadget_bin(infile,npart,ngas,noffset,mode)
    implicit none

    character(kind=1,len=*) :: infile     ! output file in GADGET binary format
    integer(kind=4) :: npart,ngas          ! number of particles, gas particles
    integer(kind=4) :: noffset             ! offset
    integer(kind=4) :: npart_in_file,ngas_in_file  ! number of particles, gas particles    
    integer(kind=4) :: mode         ! mode=1 for SnapFormat=1, mode=2 for SnapFormat=2

    integer(kind=4) :: i,j,k,l          ! Dummy variables

    integer(kind=4) :: mass_block       ! Check to see if mass block is present    
    integer(kind=4) :: read_masses      ! Check to see if mass block is present    
    character(kind=1,len=4) :: check_block

    real(kind=4), allocatable :: dummy(:)
    integer(kind=4), allocatable :: idummy(:)        

    integer(kind=4) :: ioerr
    
    read_masses=1
    mass_block=0

    write(*,*) 'Reading GADGET binary file from ',trim(infile)

    open(1,file=infile,status='old',form='unformatted')
    read(1) check_block
    if(check_block.eq.'HEAD') mode=2
    rewind 1
    
    if(mode.eq.2) read(1)
    read(1) np,massarr,time,time,flagsfr,flagfeedback,nall,&
         &     flagcooling,NumFiles,BoxSize,Omega0,OmegaLambda,HubbleParam,&
         &     flagstellarage,flagmetals,highword,flagentropy,flagdblprc,&
         &     unused

    redshift=time
    expansion=1./(1.+redshift)
    time=expansion
    
    npart_in_file=0
    
    do i=1,ntype
       npart_in_file=npart_in_file+nall(i)
       if(nall(i).gt.0.and.massarr(i).gt.0.0) read_masses=0
    end do

    ngas_in_file=nall(1)

    j=noffset

    do k=1,ntype
       do l=1,nall(k)
          j=j+1
          ptype(j)=k
       end do
    end do
    
    allocate(dummy(3*npart_in_file))
    
    if(mode.eq.2) read(1) 
    write(*,*) 'Reading positions...'
    read(1) (dummy(i),i=1,3*npart_in_file)

    j=noffset
    
    do k=1,3*npart_in_file,3
       j=j+1
       pos(1,j)=dummy(k)
       pos(2,j)=dummy(k+1)
       pos(3,j)=dummy(k+2)             
    end do

    if(mode.eq.2) read(1) 
    write(*,*) 'Reading velocities...'
    read(1) (dummy(i),i=1,3*npart_in_file)    

    j=noffset
    
    do k=1,3*npart_in_file,3
       j=j+1
       vel(1,j)=dummy(k)
       vel(2,j)=dummy(k+1)
       vel(3,j)=dummy(k+2)             
    end do

    deallocate(dummy)    

    allocate(idummy(npart_in_file))
    if(mode.eq.2) read(1) 
    write(*,*) 'Reading identities...'
    read(1) (idummy(i),i=1,npart_in_file)        

    j=noffset
    do k=1,npart_in_file
       j=j+1       
       id(j)=idummy(k)
    end do

    deallocate(idummy)

    allocate(dummy(npart_in_file))
    
    if(mode.eq.2) then
       read(1) check_block
       if(check_block.eq.'MASS') mass_block=1
    end if

    if(read_masses.eq.1.or.mass_block.eq.1) then
       write(*,*) 'Reading masses...'
       read(1) (dummy(i),i=1,npart_in_file)
       j=noffset
       do k=1,npart_in_file
          j=j+1
          mass(j)=dummy(k)
       end do
    end if

    if(read_masses.eq.0) then
       j=noffset
       write(*,*) 'Assigning masses from massarr...'
       do k=1,ntype
          do l=1,nall(k)
             j=j+1
             mass(j)=massarr(k)
          end do
       end do
    end if
    
    if(ngas.gt.0) then
       if(mode.eq.2.and.mass_block.eq.1) read(1)
       write(*,*) 'Reading internal energies...'
       read(1) (dummy(i),i=1,ngas_in_file)        
       j=noffset
       do k=1,ngas_in_file
          j=j+1
          ug(j)=dummy(k)
       end do
       if(mode.eq.2.and.mass_block.eq.1) read(1)
       write(*,*) 'Reading smoothing_lengths...'
       read(1,iostat=ioerr) (dummy(i),i=1,ngas_in_file)
       if(ioerr<0) then
          write(*,*) 'Skipping reading smoothering lengths'
          j=noffset
          do k=1,ngas_in_file
             j=j+1
             hsml(j)=dummy(k)
          end do
       endif
    end if
    deallocate(dummy)
    close(1)
    write(*,*) 'Finished reading input'

    noffset=npart_in_file
    
    return
  end subroutine read_gadget_bin

#ifdef HDF5  
  subroutine read_hdf5_header(infile,npart,ngas,mode)
    implicit none
    
    character(kind=1,len=*) :: infile  ! input file
    integer(kind=4) :: npart,ngas      ! number of particles, gas particles
    integer(kind=4) :: mode   ! mode=3 for GADGET, mode=4 for SWIFT

    integer(hsize_t), dimension(2) :: data_dims
    integer(kind=4) :: rank    

    integer(kind=4) :: i,j,k,l,m,n    ! Dummy variables
    
    integer(hid_t) :: file_id         ! File identifier
    integer(hid_t) :: group_id        ! Dataset identifier
    
    integer(kind=4) :: io_error
    character(kind=1,len=20) :: name

    integer(kind=4) :: nmembers

    integer(kind=4) :: type
    
    call h5open_f(io_error)
    
    call h5fopen_f(infile,H5F_ACC_RDONLY_F,file_id,io_error)

    call h5gopen_f(file_id,"/Header",group_id,io_error)
    
    call read_header_hdf5(group_id,mode)
    
    call h5gclose_f(group_id,io_error)

    if(mode.eq.4) then
       call h5gopen_f(file_id,"/Units",group_id,io_error)
       
       call read_units_hdf5(group_id)
       
       call h5gclose_f(group_id,io_error)
    end if
       
    call h5fclose_f(file_id,io_error)

    call h5close_f(io_error)

    npart=0

    do i=1,ntype
       npart=npart+nall(i)
    end do
    
    ngas=np(1)

    return

  end subroutine read_hdf5_header

  subroutine write_hdf5(outfile,npart,ngas,mode)
    implicit none

    character(kind=1,len=*) :: outfile     ! output file in GADGET binary format
    integer(kind=4) :: npart,ngas          ! number of particles, gas particles
    integer(kind=4) :: mode         ! mode=1 for SnapFormat=1, mode=2 for SnapFormat=2

    integer(hsize_t), dimension(2) :: data_dims
    integer(kind=4) :: rank    

    integer(kind=4) :: write_masses        ! write to block if 1    
    integer(kind=4) :: i,j,k,l,m,n    ! Dummy variables
    
    integer(hid_t) :: file_id         ! File identifier
    integer(hid_t) :: group_id        ! Dataset identifier
    
    integer(kind=4) :: io_error
    character(kind=1,len=20) :: name

    j=0
    k=0

    print *, 'np', np
    do i=1,ntype
       if(np(i).gt.0) then
          k=k+1
          if(massarr(i).gt.0.0) j=j+1
       end if
    end do

    write_masses=1

    if(j.eq.k) then
       write_masses=0
    else
       j=0
       do i=1,ntype
          if(np(i).gt.0.and.massarr(i).gt.0.0) then
             do k=j+1,j+np(i)
                mass(k)=massarr(i)
             end do
             massarr(i)=0.0
          end if
          j=np(i)
       end do       
    end if
    
    if(ismassblock) write_masses=1

    if(mode.eq.4) write_masses=1

    if(mode.eq.3) then
       write(*,*) 'Writing GADGET HDF5 file to ',trim(outfile)
    else
       write(*,*) 'Writing SWIFT HDF5 file to ',trim(outfile)
    end if
    
    call h5open_f(io_error)
    
    call h5fcreate_f(outfile,H5F_ACC_TRUNC_F,file_id,io_error)
    
    call h5gcreate_f(file_id, "/Header", group_id, io_error);
    
    call write_header_hdf5(group_id,mode)    
    
    call h5gclose_f(group_id,io_error)

    if(mode.eq.4) then
       call h5gcreate_f(file_id, "/Units", group_id, io_error);
    
       call write_units_hdf5(group_id)    
       
       call h5gclose_f(group_id,io_error)
    end if

    k=1
    do i=1,ntype
       print *, i, nall(i)
       l=k+nall(i)-1
       if(nall(i).gt.0) then
          write(name,'(a9,i1)') '/PartType',i-1
          call h5gcreate_f(file_id,name,group_id,io_error)
          rank=2
          data_dims(1)=3
          data_dims(2)=nall(i)
          name="Coordinates"
          call write_block_hdf5(group_id,1,name,rank,data_dims,3*nall(i),pos(:,k:l))
          name="Velocities"
          call write_block_hdf5(group_id,1,name,rank,data_dims,3*nall(i),vel(:,k:l))
          rank=1
          data_dims(1)=nall(i)

          name="ParticleIDs"
          call write_block_hdf5(group_id,0,name,rank,data_dims,nall(i),real(id(k:l)))
          if(write_masses.eq.1) then
             name="Masses"
             call write_block_hdf5(group_id,1,name,rank,data_dims,nall(i),mass(k:l))
          end if
          if(i.eq.1) then
             name="InternalEnergy"
             call write_block_hdf5(group_id,1,name,rank,data_dims,nall(i),ug)
             if(mode.eq.4) then
                name="SmoothingLength"
                call write_block_hdf5(group_id,1,name,rank,data_dims,nall(i),hsml)
             end if
          end if
          
          call h5gclose_f(group_id,io_error)
       end if
       k=k+nall(i)
    end do
    
    call h5fclose_f(file_id,io_error)

    call h5close_f(io_error)

    write(*,*) 'Finished writing output'

    return

  end subroutine write_hdf5

  subroutine read_hdf5(infile,npart,ngas,noffset,mode)
    implicit none

    integer(kind=4) :: mode
    integer(kind=4) :: npart,ngas,ndata
    integer(kind=4) :: noffset            ! offset
    character(kind=1,len=*) :: infile     ! output file in GADGET binary format

    integer(hsize_t), dimension(2) :: data_dims
    integer(kind=4) :: rank    

    integer(kind=4) :: i,j,k,l            ! Dummy variables
    
    integer(hid_t) :: file_id             ! File identifier
    integer(hid_t) :: group_id            ! Dataset identifier
    
    integer(kind=4) :: io_error
    character(kind=1,len=20) :: name

    real(kind=4), allocatable :: dummy(:)

    integer(kind=4) :: read_masses=0

    call h5open_f(io_error)
    
    call h5fopen_f(infile,H5F_ACC_RDONLY_F,file_id,io_error)

    call h5gopen_f(file_id,"/Header",group_id,io_error)

    call read_header_hdf5(group_id,mode)
    
    call h5gclose_f(group_id,io_error)

    if(mode.eq.4) then
       call h5gopen_f(file_id,"/Units",group_id,io_error)
       
       call read_units_hdf5(group_id)
       
       call h5gclose_f(group_id,io_error)
    end if
    
    do i=1,ntype

       j=noffset
       do l=1,np(i)
          j=j+1
          ptype(j)=i
       end do

       read_masses=1

       if(massarr(i).gt.0) read_masses=0
       
       if(np(i).gt.0) then
          
          write(name,'(a9,i1)') '/PartType',i-1
          call h5gopen_f(file_id,name,group_id,io_error)
          rank=2
          data_dims(1)=3
          data_dims(2)=nall(i)

          ndata=data_dims(1)*data_dims(2)

          allocate(dummy(ndata))

          name="Coordinates"
          call read_block_hdf5(group_id,1,name,rank,data_dims,ndata,dummy,io_error)

          j=noffset
          
          do k=1,ndata,3
             j=j+1
             pos(1,j)=dummy(k)
             pos(2,j)=dummy(k+1)
             pos(3,j)=dummy(k+2)
          end do

          name="Velocities"
          call read_block_hdf5(group_id,1,name,rank,data_dims,ndata,dummy,io_error)

          j=noffset
          
          do k=1,ndata,3
             j=j+1
             vel(1,j)=dummy(k)
             vel(2,j)=dummy(k+1)
             vel(3,j)=dummy(k+2)             
          end do
          
          deallocate(dummy)
          
          rank=1
          data_dims(1)=nall(i)
          data_dims(2)=1

          ndata=data_dims(1)*data_dims(2)
          
          allocate(dummy(ndata))         

          name="ParticleIDs"
          call read_block_hdf5(group_id,0,name,rank,data_dims,ndata,dummy,io_error)

          j=noffset
          
          do k=1,ndata
             j=j+1
             id(j)=real(dummy(k))
          end do

          if(read_masses.eq.1) then
             name="Masses"
             call read_block_hdf5(group_id,1,name,rank,data_dims,ndata,dummy,io_error)
             
             j=noffset
             
             do k=1,ndata
                j=j+1
                mass(j)=dummy(k)
             end do
          else
             j=noffset
             
             do k=1,nall(i)
                j=j+1
                mass(j)=massarr(i)
             end do
          end if
             
          if(i.eq.1) then
             name="InternalEnergy"
             call read_block_hdf5(group_id,1,name,rank,data_dims,ndata,dummy,io_error)
             j=noffset
             do k=1,ndata
                j=j+1
                ug(j)=dummy(k)
             end do
             if(mode.eq.4) then
                name="SmoothingLength"
                call read_block_hdf5(group_id,1,name,rank,data_dims,ndata,dummy,io_error)
                j=noffset
                do k=1,ndata
                   j=j+1
                   hsml(j)=dummy(k)
                end do
             end if
          end if

          deallocate(dummy)
          
          call h5gclose_f(group_id,io_error)
       end if
       noffset=j
    end do
    
    call h5fclose_f(file_id,io_error)

    call h5close_f(io_error)

    return

  end subroutine read_hdf5
    
  subroutine write_header_hdf5(handle,mode)
    implicit none

    integer(kind=4) :: mode ! 3 for GADGET, 4 for SWIFT
    integer(hsize_t), dimension(1) :: adim,data_dims
    integer(kind=4) :: rank
    integer(kind=4) :: io_error    

    integer(hid_t) :: hdf5_dataspace
    integer(hid_t) :: hdf5_attribute
    integer(hid_t) :: handle

    rank=1
    adim(1)=ntype
    data_dims(1)=ntype

    ! NumPart_ThisFile
    call h5screate_f(H5S_SIMPLE_F,hdf5_dataspace,io_error)
    call h5sset_extent_simple_f(hdf5_dataspace,rank,adim,adim,io_error)
    call h5acreate_f(handle,"NumPart_ThisFile",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,np,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    ! NumPart_Total    
    call h5screate_f(H5S_SIMPLE_F,hdf5_dataspace,io_error)
    call h5sset_extent_simple_f(hdf5_dataspace,rank,adim,adim,io_error)
    call h5acreate_f(handle,"NumPart_Total",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,nall,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    ! NumPart_Total_Highword    
    call h5screate_f(H5S_SIMPLE_F,hdf5_dataspace,io_error)
    call h5sset_extent_simple_f(hdf5_dataspace,rank,adim,adim,io_error)
    call h5acreate_f(handle,"NumPart_Total_HighWord",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,highword,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    ! MassTable
    call h5screate_f(H5S_SIMPLE_F,hdf5_dataspace,io_error)
    call h5sset_extent_simple_f(hdf5_dataspace,rank,adim,adim,io_error)
    call h5acreate_f(handle,"MassTable",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,massarr,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    adim(1)=1
    data_dims(1)=1
    
    ! Basic parameters
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Time",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,time,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    if(mode.eq.3) then
       call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
       call h5acreate_f(handle,"BoxSize",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
            & io_error)
       call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,BoxSize,data_dims,io_error)
       call h5aclose_f(hdf5_attribute,io_error)
       call h5sclose_f(hdf5_dataspace,io_error)
    else
       adim(1)=3
       data_dims(1)=3
       call h5screate_f(H5S_SIMPLE_F,hdf5_dataspace,io_error)
       call h5sset_extent_simple_f(hdf5_dataspace,rank,adim,adim,io_error)
       call h5acreate_f(handle,"BoxSize",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
            & io_error)
       call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,(/BoxSize,BoxSize,BoxSize/),data_dims,io_error)
       call h5aclose_f(hdf5_attribute,io_error)
       call h5sclose_f(hdf5_dataspace,io_error)
    end if

    adim(1)=1
    data_dims(1)=1
    
    ! Cosmological parameters
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Redshift",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,redshift,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    if(mode.eq.4) then
       call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
       call h5acreate_f(handle,"Scale-factor",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
            & io_error)
       call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,1./(1.+redshift),data_dims,io_error)
       call h5aclose_f(hdf5_attribute,io_error)
       call h5sclose_f(hdf5_dataspace,io_error)
    end if
    
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    if(mode.eq.3) then
       call h5acreate_f(handle,"Omega0",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
            & io_error)
    else
       call h5acreate_f(handle,"Omega_m",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
            & io_error)
    end if
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,Omega0,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    if(mode.eq.3) then
       call h5acreate_f(handle,"OmegaLambda",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
            & io_error)
    else
       call h5acreate_f(handle,"Omega_l",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
            & io_error)
    end if
       
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,OmegaLambda,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"HubbleParam",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,HubbleParam,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    ! Number of files per snapshot
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"NumFilesPerSnapshot",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,NumFiles,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    ! Various flags
    
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Flag_Cooling",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagcooling,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Flag_DoublePrecision",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagdp,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)
    
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Flag_Entropy_ICs",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagentropy,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)    

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Flag_Feedback",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagfeedback,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)
    
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Flag_Metals",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagmetals,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)    

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Flag_Sfr",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagsfr,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Flag_StellarAge",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagstellarage,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)
    
    if(mode.eq.4) then
       call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
       call h5acreate_f(handle,"Dimension",H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_attribute,&
            & io_error)
       call h5awrite_f(hdf5_attribute,H5T_NATIVE_INTEGER,dims,data_dims,io_error)
       call h5aclose_f(hdf5_attribute,io_error)
       call h5sclose_f(hdf5_dataspace,io_error)
    end if
    
  end subroutine write_header_hdf5

  subroutine read_header_hdf5(handle,mode)
    implicit none

    integer(kind=4) :: mode   ! 3 for GADGET, 4 for SWIFT
    integer(hsize_t), dimension(1) :: adim,data_dims
    integer(kind=4) :: rank
    integer(kind=4) :: io_error    

    integer(hid_t) :: hdf5_dataspace
    integer(hid_t) :: hdf5_attribute
    integer(hid_t) :: hdf5_link_access_id
    integer(hid_t) :: handle

    logical :: link_exist

    real(kind=4), dimension(3) :: dummy
    
    rank=1
    data_dims(1)=ntype

    call h5aopen_name_f(handle,"NumPart_ThisFile",hdf5_attribute,io_error)
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,np,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"NumPart_Total",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,nall,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"NumPart_Total_HighWord",hdf5_attribute,io_error)        
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,highword,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"MassTable",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,massarr,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    data_dims(1)=1
    
    ! Basic parameters
    call h5aopen_name_f(handle,"Time",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,time,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    if(mode.eq.3) then
       call h5aopen_name_f(handle,"BoxSize",hdf5_attribute,io_error)        
       call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,BoxSize,data_dims,io_error)
       call h5aclose_f(hdf5_attribute,io_error)
    else
       data_dims(1)=3
       call h5aopen_name_f(handle,"BoxSize",hdf5_attribute,io_error)        
       call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,dummy,data_dims,io_error)
       call h5aclose_f(hdf5_attribute,io_error)
       BoxSize=dummy(1)
    end if

    data_dims(1)=1
       
    ! Cosmological parameters
    call h5aopen_name_f(handle,"Redshift",hdf5_attribute,io_error)        
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,redshift,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Omega0",hdf5_attribute,io_error)            
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,Omega0,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"OmegaLambda",hdf5_attribute,io_error)            
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,OmegaLambda,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"HubbleParam",hdf5_attribute,io_error)            
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,HubbleParam,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    ! Number of files per snapshot
    call h5aopen_name_f(handle,"NumFilesPerSnapshot",hdf5_attribute,io_error)            
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,NumFiles,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    ! Various flags
    
    call h5aopen_name_f(handle,"Flag_Cooling",hdf5_attribute,io_error)                
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagcooling,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    if(mode.ne.3) then
       call h5aopen_name_f(handle,"Flag_Entropy_ICs",hdf5_attribute,io_error)    
       call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagentropy,data_dims,io_error)
       call h5aclose_f(hdf5_attribute,io_error)
    end if
    
    call h5aopen_name_f(handle,"Flag_Feedback",hdf5_attribute,io_error)                
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagfeedback,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Flag_Metals",hdf5_attribute,io_error)                    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagmetals,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Flag_Sfr",hdf5_attribute,io_error)                
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagsfr,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Flag_StellarAge",hdf5_attribute,io_error)                    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_INTEGER,flagstellarage,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    return
    
  end subroutine read_header_hdf5

  subroutine read_units_hdf5(handle)
    implicit none

    integer(hsize_t), dimension(1) :: data_dims
    integer(kind=4) :: rank
    integer(kind=4) :: io_error    

    integer(hid_t) :: hdf5_attribute
    integer(hid_t) :: handle

    rank=1
    data_dims(1)=1

    call h5aopen_name_f(handle,"Unit length in cgs (U_L)",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_length,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Unit mass in cgs (U_M)",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_mass,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Unit time in cgs (U_t)",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_time,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Unit current in cgs (U_I)",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_current,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)

    call h5aopen_name_f(handle,"Unit temperature in cgs (U_T)",hdf5_attribute,io_error)    
    call h5aread_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_current,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    
    return
    
  end subroutine read_units_hdf5

  subroutine write_units_hdf5(handle)
    implicit none

    integer(kind=4) :: mode ! 3 for GADGET, 4 for SWIFT
    integer(hsize_t), dimension(1) :: adim,data_dims
    integer(kind=4) :: rank
    integer(kind=4) :: io_error    

    integer(hid_t) :: hdf5_dataspace
    integer(hid_t) :: hdf5_attribute
    integer(hid_t) :: handle

    rank=1
    adim(1)=1
    data_dims(1)=1
    
    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Unit length in cgs (U_L)",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_length,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Unit mass in cgs (U_M)",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_mass,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Unit time in cgs (U_t)",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_time,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Unit current in cgs (U_I)",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_current,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)

    call h5screate_f(H5S_SCALAR_F,hdf5_dataspace,io_error)
    call h5acreate_f(handle,"Unit temperature in cgs (U_T)",H5T_NATIVE_DOUBLE,hdf5_dataspace,hdf5_attribute,&
         & io_error)
    call h5awrite_f(hdf5_attribute,H5T_NATIVE_DOUBLE,unit_temperature,data_dims,io_error)
    call h5aclose_f(hdf5_attribute,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)
    
    return
    
  end subroutine write_units_hdf5
  
  subroutine write_block_hdf5(handle,itype,tag,rank,data_dims,ndata,array)
    implicit none

    integer(kind=4) :: itype
    integer(hsize_t), dimension(2) :: data_dims
    integer(kind=4) :: rank
    integer(kind=4) :: io_error    

    integer(hid_t) :: hdf5_dataspace
    integer(hid_t) :: hdf5_dataset
    integer(hid_t) :: handle

    character(kind=1,len=*) :: tag

    integer(kind=4) :: ndata
    real(kind=4), dimension(ndata) :: array

    call h5screate_simple_f(rank,data_dims,hdf5_dataspace,io_error)
    if(itype.eq.1) then
       call h5dcreate_f(handle,tag,H5T_NATIVE_REAL,hdf5_dataspace,hdf5_dataset,io_error)
       call h5dwrite_f(hdf5_dataset,H5T_NATIVE_REAL,array,data_dims,io_error)
    else
       call h5dcreate_f(handle,tag,H5T_NATIVE_INTEGER,hdf5_dataspace,hdf5_dataset,io_error)
       call h5dwrite_f(hdf5_dataset,H5T_NATIVE_INTEGER,int(array),data_dims,io_error)       
    end if
    call h5dclose_f(hdf5_dataset,io_error)
    call h5sclose_f(hdf5_dataspace,io_error)
    
  end subroutine write_block_hdf5

  subroutine read_block_hdf5(handle,itype,tag,rank,data_dims,ndata,array,io_error)
    implicit none

    integer(kind=4) :: itype
    integer(hsize_t), dimension(2) :: data_dims
    integer(kind=4) :: rank
    integer(kind=4) :: io_error    

    integer(hid_t) :: hdf5_dataset
    integer(hid_t) :: handle

    character(kind=1,len=*) :: tag

    integer(kind=4) :: ndata
    real(kind=4), dimension(ndata) :: array
    integer(kind=4), dimension(ndata) :: iarray    

    call h5dopen_f(handle,tag,hdf5_dataset,io_error)
    if(io_error.eq.-1) return
    if(itype.eq.1) then
       call h5dread_f(hdf5_dataset,H5T_NATIVE_REAL,array,data_dims,io_error)
    else
       call h5dread_f(hdf5_dataset,H5T_NATIVE_INTEGER,iarray,data_dims,io_error)
       array=real(iarray)
    end if
    call h5dclose_f(hdf5_dataset,io_error)
    
  end subroutine read_block_hdf5
#endif

  subroutine sort_particles_by_ptype(npart,ngas,npall)
    implicit none
    
    integer(kind=4), dimension(ntype) :: npall
    integer(kind=4) :: npart,ngas
    integer(kind=4), dimension(ntype) :: noffset,n_in_type
    integer(kind=4) :: i,j,k
    
    real(kind=4), allocatable :: dummy(:,:)

    j=0
    
    do i=1,ntype
       noffset(i)=j
       j=j+npall(i)
    end do

    allocate(dummy(3,npart))
    
    ! Positions
    do i=1,ntype
       n_in_type(i)=0
    end do
    
    do i=1,npart
       n_in_type(ptype(i))=n_in_type(ptype(i))+1
       j=n_in_type(ptype(i))+noffset(ptype(i))
       dummy(1,j)=pos(1,i)
       dummy(2,j)=pos(2,i)
       dummy(3,j)=pos(3,i)     
    end do

    do i=1,npart
       do j=1,3
          pos(j,i)=dummy(j,i)
       end do
    end do
    
    ! Velocities
    do i=1,ntype
       n_in_type(i)=0
    end do
    
    do i=1,npart
       n_in_type(ptype(i))=n_in_type(ptype(i))+1
       j=n_in_type(ptype(i))+noffset(ptype(i))
       dummy(1,j)=vel(1,i)
       dummy(2,j)=vel(2,i)
       dummy(3,j)=vel(3,i)     
    end do
    
    do i=1,npart
       do j=1,3
          vel(j,i)=dummy(j,i)
       end do
    end do
    
    ! Masses
    do i=1,ntype
       n_in_type(i)=0
    end do
    
    do i=1,npart
       n_in_type(ptype(i))=n_in_type(ptype(i))+1
       j=n_in_type(ptype(i))+noffset(ptype(i))
       dummy(1,j)=mass(i)
    end do

    do i=1,npart
       mass(i)=dummy(1,i)
    end do
    
    ! Internal Energies
    do i=1,ntype
       n_in_type(i)=0
    end do
    
    do i=1,npart
       if(ptype(i).eq.1) then
          n_in_type(ptype(i))=n_in_type(ptype(i))+1
          j=n_in_type(ptype(i))+noffset(ptype(i))
          dummy(1,j)=ug(i)
       end if
    end do
    
    do i=1,n_in_type(1)
       ug(i)=dummy(1,i)
    end do
    
    deallocate(dummy)

    ! Overwrite IDs - avoids issues with e.g. GalIC outputs

    do i=1,npart
       id(i)=i
    end do

    do i=1,ntype
       np(i)=npall(i)
       nall(i)=npall(i)
    end do
    
    return
    
  end subroutine sort_particles_by_ptype

  subroutine sample_particles_by_id(npart,npall,indx)
    implicit none

    integer(kind=4) :: npart
    integer(kind=4), dimension(ntype) :: npall
    integer(kind=4), dimension(npart) :: indx,jndx

    integer(kind=4) :: i,j,k
    
    real(kind=4), allocatable :: dummy(:,:)

    allocate(dummy(3,npart))

    do i=1,npart
       dummy(1,i)=indx(i)
    end do
    
    call indexx(npart,dummy(1,:),jndx)
    
    do i=1,npart
       j=jndx(i)
       pos(1,i)=pos(1,id(indx(j)))
       pos(2,i)=pos(2,id(indx(j)))
       pos(3,i)=pos(3,id(indx(j)))
       vel(1,i)=vel(1,id(indx(j)))
       vel(2,i)=vel(2,id(indx(j)))
       vel(3,i)=vel(3,id(indx(j)))
       mass(i)=mass(id(indx(j)))
       ptype(i)=ptype(id(indx(j)))
       if(ptype(i).eq.1) ug(i)=ug(id(indx(j)))
       id(i)=i
    end do
    
    deallocate(dummy)

    do i=1,ntype
       np(i)=0
    end do
    
    do i=1,npart
       np(ptype(i))=np(ptype(i))+1
    end do

    j=1
    do i=1,ntype
       write(*,*) 'Species',i,' :',nall(i),np(i)
       npall(i)=np(i)
       nall(i)=npall(i)
    end do

    return
    
  end subroutine sample_particles_by_id
  

end module io



