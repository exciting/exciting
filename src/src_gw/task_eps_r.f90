      
subroutine task_eps_r

    use modinput
    use modmain
    use modgw
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
    use mod_rpath
            
    implicit none
    integer(4) :: ikp, iq, fid, ik
    real(8)    :: t0, t1
    integer(4) :: recl
    integer(4) :: im, iom, npt, ir, ir0, theta0
    complex(8) :: vc
    
    complex(8), allocatable :: wfmb(:,:)
    complex(8), allocatable :: tvec(:), tmat(:,:)
    complex(8), allocatable :: eps_q (:,:), eps_r(:)

    complex(8), external :: zdotc

    !===========================================================================
    ! Initialization
    !===========================================================================

    ! Frequency point
    iom = 1

    ! Position of the first real space point
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
    
    if (.not.input%gw%rpmat) then
      !========================================================
      ! calculate momentum matrix elements and store to a file
      !========================================================
      call calcpmatgw
    end if
    
    ! occupancy dependent BZ integration weights
    call kintw

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
        write(*,*) "ERROR(task_eps_r): Unknown path type!", trim(input%gw%rpath)
        stop

    end select

    npt = rpath%nptot

    allocate(wfmb(matsizmax,npt))
    allocate(tvec(matsizmax))
    allocate(eps_r(npt))
    eps_r(:) = 0.d0
    
    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================    

#ifdef MPI
    call set_mpi_group(kqset%nkpt)
    call mpi_set_range(nproc_row, &
    &                  myrank_row, &
    &                  kset%nkpt, 1, &
    &                  iqstart, iqend)
#else
    iqstart = 1
    iqend = kqset%nkpt
#endif
    iomstart = 1
    iomend = freq%nomeg
    
    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if
    
    ! each process does a subset
    do iq = iqstart, iqend
    
      write(*,*)
      write(*,*) '(task_eps_r): q-point cycle, iq = ', iq
    
      Gamma = gammapoint(kqset%vqc(:,iq))
          
      !========================================
      ! Calculate interstitial basis functions
      !========================================
      matsiz = locmatsiz+Gqset%ngk(1,iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
    
      !======================================
      ! Calculate the bare Coulomb potential
      !======================================
      call calcbarcmb(iq)

      !========================================
      ! Set v-diagonal MB and reduce its size
      !========================================
      call setbarcev(input%gw%barecoul%barcevtol)
      
      !===================================
      ! Calculate the dielectric function
      !===================================
      call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)

      ! dielectric matrix      
      call calcepsilon(iq,iomstart,iomend)

      ! e = (1-vP) - 1
      do im = 1, mbsiz
        epsilon(im,im,iomstart:iomend) = epsilon(im,im,iomstart:iomend)-zone
      end do

if (.false.) then
      ! inverse dielectic matrix
      call calcinveps(iomstart,iomend)
      do im = 1, mbsiz
        epsilon(im,im,iom) = epsilon(im,im,iom)+zone
      end do
end if

      allocate(tmat(mbsiz,matsiz))
      call zgemm('n','c', &
      &          mbsiz,matsiz,mbsiz, &
      &          zone, &
      &          epsilon(:,:,iom),mbsiz, &
      &          vmat,matsiz, &
      &          zzero,tmat,mbsiz)

      allocate(eps_q(matsiz,matsiz))
      call zgemm('n','n', &
      &          matsiz,matsiz,mbsiz, &
      &          zone, &
      &          vmat,matsiz, &
      &          tmat,mbsiz, &
      &          zzero,eps_q,matsiz)
      deallocate(tmat)

      !===================================
      ! Real-space representation
      !===================================
      select case (trim(input%gw%rpath))
        case("atoms")
          call calc_mb_functions(iq,npt,wfmb)
          call zgemv('n',mbsiz,mbsiz,zone,eps_q,mbsiz, &
          &          wfmb(:,ir0),1,zzero,tvec,1)

        case("rad")
          ! radial direction
          call calc_radial_wfmb(iq,npt,wfmb)
          call zgemv('n',mbsiz,mbsiz,zone,eps_q,mbsiz, &
          &          wfmb(:,ir0),1,zzero,tvec,1)

        case("azi")
          ! azimuthal direction
          call calc_azimuthal_wfmb(iq,ir0,npt,wfmb)
          call zgemv('n',mbsiz,mbsiz,zone,eps_q,mbsiz, &
          &          wfmb(:,theta0),1,zzero,tvec,1)

      end select
      deallocate(eps_q)

      do ir = 1, npt
        eps_r(ir) = eps_r(ir)+zdotc(mbsiz,wfmb(:,ir),1,tvec,1)/dble(kqset%nkpt)
      end do

      call delete_coulomb_potential
      call delete_dielectric_function(Gamma)
      if (allocated(kcw)) deallocate(kcw)
      if (allocated(unw)) deallocate(unw)
      
      ! clean unused data
      if (allocated(mpwipw)) deallocate(mpwipw)
      if (allocated(barc)) deallocate(barc)
      
    end do ! iq

    deallocate(tvec,wfmb)

#ifdef MPI
    call mpi_sum_array(0,eps_r,npt,mycomm_row)
#endif

    if (rank==0) then
      write(*,*) "ir0, r0=", ir0, spr(ir0,1)
      open(77,File="eps-r-plot.dat")
      do ir = 1, npt
        write(77,'(i8,2f18.6)') ir, eps_r(ir)
      end do
      close(77)
    end if
    deallocate(eps_r)

    !-------------------
    ! Final clean
    !-------------------
    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)
    
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

