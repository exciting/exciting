      
subroutine task_chi0_r

    use modinput
    use modmain,               only : zzero, evalsv, efermi
    use modgw
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
    use mod_rpath
            
    implicit none
    integer(4) :: ikp, iq, fid, ik
    real(8)    :: t0, t1
    integer(4) :: recl
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
    call kintw
    
    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================    

#ifdef MPI
    call set_mpi_group(kqset%nkpt)
    call mpi_set_range(nproc_row, &
    &                  myrank_row, &
    &                  kqset%nkpt, 1, &
    &                  iqstart, iqend)
#else
    iqstart = 1
    iqend = kqset%nkpt
#endif
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

    if (myrank==0) then
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
      call calcepsilon(iq,iomstart,iomend)
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
    call mpi_sum_array(0,chi0_r,npt,mycomm_row)
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
      call getunit(fid)
      open(fid,File='CHI0.OUT',form='UNFORMATTED',status='REPLACE')
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
        &               chi0(1,1,iomstart,iq),(/matsizmax,matsizmax,iomstart:iomend/))
#endif
        call getunit(fid)
        inquire(iolength=recl) chi0(:,:,:,iq)
        open(fid,File="CHI0.OUT",Action='Write',Form='Unformatted',&
        &    Access='Direct',Status='Old',Recl=recl)
        write(fid,rec=iq) chi0(:,:,:,iq)
        close(fid)
      end if ! rank
      call barrier
    end do ! iq
end if    
    
    deallocate(chi0)
    if (allocated(evalsv)) deallocate(evalsv)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return
end subroutine
