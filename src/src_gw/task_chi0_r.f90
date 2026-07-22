      
subroutine task_chi0_r
    use calculate_dielectric_function, only: calcepsilon, epsilon_indexes
    use modinput
    use mod_large_io, only: inquire_large, open_direct_unformatted_large
    use modmain,               only : zzero, efermi
    use modmpi, only: distribute_loop, mpiglobal, rank
    use modgw
    use mod_mpi_gw
    use modmpi, only: rank, mpiglobal, barrier
    use mod_hdf5
    use mod_rpath
    use mod_coulomb_potential, only: barc
    use mod_bands, only: evalfv, numin, nstdf
    use precision, only: i32, long_int
            
    implicit none
    integer(4) :: ikp, iq, ik
    integer(i32) :: fid
    real(8)    :: t0, t1
    integer(long_int) :: recl
    integer :: im, iom, npt, ir, ir0, theta0
    
    ! mapping array 
    integer, allocatable :: iq2rank(:)
    
    complex(8), allocatable :: chi0(:,:,:,:)
    complex(8), allocatable :: wfmb(:,:), tvec(:), chi0_r(:)
    
    complex(8), external :: zdotc

    !===========================================================================
    ! Initialization
    !===========================================================================

    ! Frequency point
    iom = 1

    ! radial grid point
    ir0 = input%gw%iik

    ! theta
    theta0 = input%gw%jjk

    !----------------------------------------------
    ! Store all important results to the hdf5 file
    !----------------------------------------------
#ifdef _HDF5_
    call hdf5_initialize()
    fgwh5 = "gw_output.h5"
    if (rank==0) then
      call hdf5_create_file(fgwh5)
      call hdf5_create_group(fgwh5,"/","parameters")
      if (rank==0) call write_gw_parameters_hdf5
      call hdf5_create_group(fgwh5,"/","kpoints")
    end if
#endif    
    
    ! prepare GW global data
    call init_gw
    
    ! clean not used anymore global exciting variables
    call clean_gndstate
    
    ! occupancy dependent BZ integration weights
    call kintw()
    
    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================    

    call distribute_loop( mpiglobal, kqset%nkpt, iqstart, iqend )
    iomstart = 1
    iomend = freq%nomeg

    ! iq <--> MPI rank mapping
    allocate(iq2rank(kqset%nkpt))
    iq2rank = -1
    do iq = iqstart, iqend
      iq2rank(iq) = rank
    end do
    
    if (allocated(chi0)) deallocate(chi0)
    allocate(chi0(matsizmax,matsizmax,iomstart:iomend,iqstart:iqend))
    chi0(:,:,:,:) = 0.d0
    
    ! real space definition
    select case (trim(input%gw%rpath))
      case("atoms")
        call init_rpath(rpath,input%gw%at1,input%gw%at2)

      case("rad")
        call init_radial_path(rpath,input%gw%at1)

      case("azi")
        call init_azimuthal_path(rpath,input%gw%at1,ir0)

      case default
        write(*,*)
        write(*,*) "ERROR(task_chi0_r): Unknown path type!", trim(input%gw%rpath)
        stop

    end select

    npt = rpath%nptot

    allocate(wfmb(matsizmax,npt))
    allocate(tvec(matsizmax))
    allocate(chi0_r(npt))
    chi0_r(:) = 0.d0

    if (rank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if
    
    ! each process does a subset
    do iq = iqstart, iqend
    
      write(*,*)
      write(*,*) '(task_chi0_r): q-point cycle, iq = ', iq
      
      Gamma = .false.
    
      !========================================
      ! Calculate interstitial basis functions
      !========================================
      matsiz = locmatsiz+Gqset%ngk(1,iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
    
      !======================================
      ! Coulomb potential (not used)
      !======================================
      mbsiz = matsiz
      if (allocated(barc)) deallocate(barc)
      allocate(barc(matsiz,mbsiz))
      barc(:,:) = 0.d0
      do im = 1, matsiz
        barc(im,im) = zone
      end do
            
      !===================================
      ! Calculate the chi0 function
      !===================================
      call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
      call calcepsilon(iq, epsilon_indexes( &
                      indexes_parallelization( 1, kqset%nkpt, 1, kqset%nkpt ), &
                      indexes_parallelization( numin, nstdf, numin, nstdf ), &
                      indexes_parallelization( iomstart, iomend, iomstart, iomend ) ) &
                      )
      do im = 1, mbsiz
        epsilon(im,im,iomstart:iomend) = epsilon(im,im,iomstart:iomend)-zone
      end do
    
      ! save to local array: e = 1-vP
      chi0(1:mbsiz,1:mbsiz,:,iq) = -epsilon(:,:,:)
      
      !===================================
      ! Real-space representation
      !===================================
      select case (trim(input%gw%rpath))
        case("atoms")
          call calc_mb_functions(iq,npt,wfmb)
          call zgemv('n',mbsiz,mbsiz,zone,epsilon(:,:,iom),mbsiz, &
          &          wfmb(:,ir0),1,zzero,tvec,1)
          
        case("rad")
          ! radial direction
          call calc_radial_wfmb(iq,npt,wfmb)
          call zgemv('n',mbsiz,mbsiz,zone,epsilon(:,:,iom),mbsiz, &
          &          wfmb(:,ir0),1,zzero,tvec,1)

        case("azi")
          ! azimuthal direction
          call calc_azimuthal_wfmb(iq,ir0,npt,wfmb)
          call zgemv('n',mbsiz,mbsiz,zone,epsilon(:,:,iom),mbsiz, &
          &          wfmb(:,theta0),1,zzero,tvec,1)

      end select

      do ir = 1, npt
        chi0_r(ir) = chi0_r(ir)+zdotc(mbsiz,wfmb(:,ir),1,tvec,1)/dble(kqset%nkpt)
      end do
      
      call delete_dielectric_function(Gamma)
      if (allocated(kcw)) deallocate(kcw)
      if (allocated(unw)) deallocate(unw)
      
      ! clean unused data
      if (allocated(mpwipw)) deallocate(mpwipw)
      if (allocated(barc)) deallocate(barc)
      
    end do ! iq
    
    deallocate(tvec,wfmb)
    
#ifdef MPI
    call mpi_sum_array( chi0_r, mpiglobal, .false. )
#endif

    if (rank==0) then
      write(*,*) "ir0, r0=", ir0, spr(ir0,1)
      open(77,File="chi0-r-plot.dat")
      do ir = 1, npt
        write(77,'(i8,2f18.6)') ir, chi0_r(ir)
      end do
      close(77)
    end if
    deallocate(chi0_r)
    
    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)

if (.false.) then    
    !===================================
    ! Write \chi_0 to file
    !===================================
    
    ! overwrite existing files
    if (rank==0) then
      open(newunit=fid,File='CHI0.OUT',form='UNFORMATTED',status='REPLACE')
      close(fid)
    endif
    call barrier
    
    do iq = 1, kqset%nkpt
      if (rank==iq2rank(iq)) then
#ifdef _HDF5_
        write(cik,'(I4.4)') ik
        path = "/qpoints/"//trim(adjustl(cik))
        if (.not.hdf5_exist_group(fgwh5,"/qpoints",cik)) &
        &  call hdf5_create_group(fgwh5,"/qpoints",cik)
        call hdf5_write(fgwh5,path,"chi0", &
        &               chi0(1,1,iomstart,iq),(/matsizmax,matsizmax,iomend-iomstart+1/))
#endif
        call inquire_large( recl, chi0(:,:,:,iq) )
        call open_direct_unformatted_large( fid, "CHI0.OUT", "write", recl, "old" )
        write(fid,rec=iq) chi0(:,:,:,iq)
        close(fid)
      end if ! rank
      call barrier
    end do ! iq
end if    
    
    deallocate(chi0)
    if (allocated(evalfv)) deallocate(evalfv)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return
end subroutine
