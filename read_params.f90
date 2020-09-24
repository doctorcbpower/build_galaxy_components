module read_params
  use structure
contains
  subroutine set_parameters(infile,npart,outfile,glass_file,istabulated,tabfile,ispoisson,mode,mode_io)
    implicit none
    
    integer(kind=4) :: mode                        ! mode=0,1 for disc, halo
    integer(kind=4) :: mode_io                     ! mode_io=0,1,2 for SnapFormat=1,2, HDF5
    character(kind=1,len=*) :: infile              ! input parameter file
    character(kind=1,len=*) :: glass_file          ! glass cube
    character(kind=1,len=*) :: outfile             ! output ICs
    character(kind=1,len=*) :: tabfile             ! tabulated profile
    character(kind=1,len=132) :: disc_file,halo_file,nbody_file  ! Name of disc and halo output files
    integer(kind=4) :: npart,ndisc,nhalo           ! Number of particles,in disc and halo
    integer(kind=4) :: snap_format                 ! GADGET binary format
    logical :: istabulated                         ! Use tabulated data?
    logical :: ispoisson                           ! Use poisson distribution?
    ! Namelists
    namelist /halo/ m200,v200,c200,sigma200,qflat,rcore  ! Halo properties - M200, c200, V200, sigma200,flattening parameter,core radius
    namelist /disc/ mdisc,sigma0,fdisc,mstar,fhole,fdref  ! Disc mass fraction, central surface mass density,
                                                          ! disc scale length, stellar disc mass fraction,
                                                          ! disc hole radius, reference length of disc
    namelist /bulge/ mbulge,fbulge     ! Bulge mass fraction, bulge scale length fraction
    namelist /gas_disc/ temp,ndisc     ! Disc temperature, number of particles in disc
    namelist /gas_halo/ fbaryon,nhalo,lambda_b,v_radial ! Baryon fraction, 
                                       ! number of particles in halo,Bullock spin parameter,
                                       ! radial infall velocity (as fraction of V200)
    namelist /black_hole/ mbh          ! Central black hole mass fraction
    namelist /nbody/ ncomp,npart_comp,fmass_comp,form_comp,rscale_comp,rtrunc_comp,raniso_comp,istaper !
    namelist /output/ disc_file,halo_file,nbody_file,snap_format   ! Disc and halo files,
                                                        ! binary format, HDF5 on or off?
    namelist /misc/ glass_file,ispoisson,tabfile,istabulated
    
    npart=0
    
    if(mode.eq.0) then ! mode=0 is for gas disc
       open(1,file=infile,status='old')
       rewind(1)
       read(1,nml=halo)
       rewind(1)
       read(1,nml=disc)
       rewind(1)
       read(1,nml=bulge)
       rewind(1)
       read(1,nml=gas_disc)
       rewind(1)
       read(1,nml=black_hole)
       rewind(1)
       read(1,nml=output)
       rewind(1)
       read(1,nml=misc)        
       close(1)
       npart=ndisc
       outfile=trim(disc_file)
    else if (mode.eq.1) then ! mode=1 is for gas halo
        open(1,file=infile,status='old')
        rewind(1)
        read(1,nml=halo)
        rewind(1)
        read(1,nml=disc)
        rewind(1)
        read(1,nml=bulge)
        rewind(1)
        read(1,nml=gas_halo)
        rewind(1)
        read(1,nml=black_hole)
        rewind(1)
        read(1,nml=output)
        rewind(1)
        read(1,nml=misc)
        close(1)
        npart=nhalo
        outfile=trim(halo_file)
    else if (mode.eq.2) then ! mode=2 is for nbody components
        open(1,file=infile,status='old')
        rewind(1)
        read(1,nml=halo)
        rewind(1)
        read(1,nml=disc)
        rewind(1)
        read(1,nml=bulge)
        rewind(1)
        read(1,nml=black_hole)
        rewind(1)
        read(1,nml=nbody)
        rewind(1)
        read(1,nml=output)
        rewind(1)
        read(1,nml=misc)
        close(1)
        npart=sum(npart_comp)
        outfile=trim(nbody_file)
    else 
        stop 'mode not supported'    
    end if 

    !if(mode.eq.0) then        ! mode=0 is for gas disc
    !   npart=ndisc
    !   outfile=trim(disc_file)
    !else if (mode.eq.1) then  ! mode=1 is for gas halo
    !    npart=nhalo
    !    outfile=trim(halo_file)
    !else if (mode.eq.2) then  ! mode=2 is for nbody components
    !    npart=sum(npart_comp)
    !    outfile=trim(nbody_file)
    !else
    !    stop 'mode not supported'
    !end if

    if(snap_format.gt.4) snap_format=1

    mode_io=snap_format
    
    if(mode_io.eq.1) then
       outfile=trim(outfile)//'.snp1.gdt'
    else if(mode_io.eq.2) then
       outfile=trim(outfile)//'.snp2.gdt'       
    else if(mode_io.eq.3) then
       outfile=trim(outfile)//'.gdt.hdf5'
    else if(mode_io.eq.4) then
       outfile=trim(outfile)//'.swft.hdf5'       
    end if

    return
    
  end subroutine set_parameters
end module read_params
