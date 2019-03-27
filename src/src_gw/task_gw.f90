!BOP
!
!!ROUTINE: \verb"task_gw"
!
!!INTERFACE:
!
subroutine task_gw()
!
!!DESCRIPTION:
!
! This subroutine performs one GW cycle and calculates the corresponding
! quasiparticle energies.
!
!!USES:
    use modinput
    use modmain,               only: zzero, efermi
    use modgw
    use mod_coulomb_potential
    use mod_vxc,               only: vxcnn
    use mod_mpi_gw
    use m_getunit

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ikp, iq, fid, ik
    real(8)    :: t0, t1
    integer(4) :: recl
    integer    :: im
    complex(8) :: vc
    integer(4) :: Nk
    real(8)    :: omega_BZ, Vk, beta, sxdiv
    character(80) :: frmt

    integer(4) :: iom, ib
    real(8) :: w, sRe, sIm, div
    complex(8) :: dsc
    real(8), allocatable :: sf(:)

!!REVISION HISTORY:
!
! Created Nov 2013 by (DIN)
!
!EOP
!BOC
    !===========================================================================
    ! Initialization
    !===========================================================================

    ! prepare GW global data
    call init_gw()

    !=================================================
    ! Calculate the diagonal matrix elements of the
    ! DFT exchange-correlation potential
    !=================================================
    ! it is better to do it here to deallocate cfunir and vxcir arrays
    call timesec(t0)
    call calcvxcnn()
    call timesec(t1)

    ! clean not used anymore global exciting variables
    ! call clean_gndstate ! disabled since it's used later in task_second_variation

    if (input%gw%taskname.ne.'g0w0_x') then
      if (.not.input%gw%rpmat) then
        !========================================================
        ! calculate momentum matrix elements and store to a file
        !========================================================
        call calcpmatgw()
      end if
    end if

    ! occupancy dependent BZ integration weights
    call kintw()

    !---------------------------------------
    ! treatment of singularities at G+q->0
    !---------------------------------------
    select case (trim(input%gw%barecoul%cutofftype))

        case('none')

            select case (trim(input%gw%selfenergy%singularity))
              case('mpb')
                ! Auxiliary function method
                call setsingc
              case('crg')
                ! Auxiliary function method
                call calc_q0_singularities
              case('avg')
                ! Spherical average
                call vcoul_q0_3d(kqset%nkpt, singc2)
              case('rim')
                ! Spherical average
              case default
                call calc_q0_singularities
            end select

        case('0d')
            call vcoul_q0_0d(singc2)

        case('1d')
            call vcoul_q0_1d(kqset%nkpt, singc2)

        case('2d')
            call vcoul_q0_2d(kqset%nkpt, singc2)

    end select

    ! initialize self-energy arrays
    call init_selfenergy(ibgw, nbgw, kset%nkpt)

    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================

#ifdef MPI
    call set_mpi_group(kqset%nkpt)
    call mpi_set_range(nproc_row, &
    &                  myrank_row, &
    &                  kqset%nkpt, 1, &
    &                  iqstart, iqend)
    call mpi_set_range(nproc_col, &
    &                  myrank_col, &
    &                  freq%nomeg, 1, &
    &                  iomstart, iomend, &
    &                  iomcnt, iomdsp)
    ! write(*,*) "myrank_row, iqstart, iqend =", myrank_row, iqstart, iqend
    ! write(*,*) "myrank_col, iomstart, iomend =", myrank_col, iomstart, iomend
    ! write(*,*) 'iomcnt: ', iomcnt(0:nproc_col-1)
    ! write(*,*) 'iomdsp: ', iomdsp(0:nproc_col-1)
#else
    iqstart = 1
    iqend = kqset%nkpt
    iomstart = 1
    iomend = freq%nomeg
#endif

    if (myrank==0) then
      call boxmsg(fgw,'=','GW cycle')
      call flushifc(fgw)
    end if

    ! each process does a subset
    do iq = iqstart, iqend

      if (rank==0) write(fgw,*) '(task_gw): q-point cycle, iq = ', iq

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

      !===============================
      ! Calculate \Sigma^{x}_{kn}(q)
      !===============================
      call calcselfx(iq)

      if (input%gw%taskname.ne.'g0w0_x') then
        !========================================
        ! Set v-diagonal MB and reduce its size
        !========================================
        if (vccut) then
          mbsiz = matsiz
          if (allocated(barc)) deallocate(barc)
          allocate(barc(matsiz,mbsiz))
          do im = 1, matsiz
            vc = cmplx(barcev(im),0.d0,8)
            barc(:,im) = vmat(:,im)*sqrt(vc)
          end do
        else
          call setbarcev(input%gw%barecoul%barcevtol)
        end if
        call delete_coulomb_potential
        !===================================
        ! Calculate the dielectric function
        !===================================
        call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
        select case (trim(input%gw%scrcoul%scrtype))
          case('ppm','PPM')
            call calcepsilon_ppm(iq,iomstart,iomend)
          case default
            call calcepsilon(iq,iomstart,iomend)
            !==========================================
            ! Calculate the screened Coulomb potential
            !==========================================
            call calcinveps(iomstart,iomend)
        end select
        !========================================
        ! Calculate the q-dependent self-energy
        !========================================
        call calcselfc(iq)
        call delete_dielectric_function(Gamma)
        if (allocated(kcw)) deallocate(kcw)
        if (allocated(unw)) deallocate(unw)
      end if

      ! clean unused data
      if (allocated(mpwipw)) deallocate(mpwipw)
      if (allocated(barc)) deallocate(barc)

    end do ! iq

    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)

#ifdef MPI
    if ((nproc_row>1).and.(myrank_col==0)) then
      call mpi_sum_array(0,selfex,nbandsgw,kset%nkpt,mycomm_row)
      if (input%gw%taskname.ne.'g0w0_x') then
        ! G0W0 and GW0 approximations
        call mpi_sum_array(0,selfec,nbandsgw,freq_selfc%nomeg,kset%nkpt,mycomm_row)
        if (input%gw%taskname=='cohsex') then
          call mpi_sum_array(0,sigsx,nbandsgw,kset%nkpt,mycomm_row)
          call mpi_sum_array(0,sigch,nbandsgw,kset%nkpt,mycomm_row)
        end if
      end if ! selfec
    endif
#endif

    !===============================================================================
    ! output block
    !===============================================================================

    if (myrank == 0) then

      if ((input%gw%taskname /= 'g0w0_x') .and. (input%gw%selfenergy%method == "ac")) then
        ! Analytical continuation of the correlation self-energy from the complex to the real frequency axis
        call plot_selfc_iw()
        call calcselfc_ac()
      end if

 !$OMP critical

      !===============================
      ! Write self-energies to files
      !===============================
      call write_selfenergy(ibgw,nbgw,kset%nkpt,freq_selfc%nomeg)
      call plot_selfc()

      !=======================================
      ! Calculate the quasiparticle energies
      !=======================================

      ! KS band structure
      evalks(ibgw:nbgw,:) = evalfv(ibgw:nbgw,:)
      call bandstructure_analysis('KS', &
          ibgw,nbgw,kset%nkpt,evalks(ibgw:nbgw,:),efermi)

      ! solve QP equation
      call calcevalqp()
      call plot_spectral_function()
      call write_qp_energies('EVALQP.DAT')

      ! G0W0 QP band structure
      select case (input%gw%taskname)

        case('g0w0_x')
          call bandstructure_analysis('G0W0_X', &
              ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)

        case('cohsex')
          call bandstructure_analysis('COHSEX', &
              ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)

        case('g0w0','gw0')
          call bandstructure_analysis('G0W0', &
              ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)

      end select

!$OMP end critical

    end if ! myrank

    !--------------------------------------------------------
    ! Calculate quasiparticle energies in GW0 approximation
    !--------------------------------------------------------
    if (input%gw%taskname=='gw0') then

      ! self-consistent cycle
      call calcscgw0

      ! print GW0 QP band structure
      if (myrank==0) then
        call timesec(t0)
        !----------------------------------------
        ! Write quasi-particle energies to file
        !----------------------------------------
        call write_qp_energies('EVALQP-GW0.DAT')
        !call bandanalysis('GW0',ibgw,nbgw,evalqp(ibgw:nbgw,:),eferqp)
        call bandstructure_analysis('GW0', ibgw, nbgw, kset%nkpt, evalqp(ibgw:nbgw,:), eferqp)
        call timesec(t1)
        time_io = time_io+t1-t0
      end if
    end if

    if (myrank==0) then
      !----------------------------------------
      ! Save QP energies into binary file
      !----------------------------------------
      call putevalqp('EVALQP.OUT')
    end if ! myrank

    !-----------------------------------------
    ! Second-variational treatment if needed
    !-----------------------------------------
    if (associated(input%groundstate%spin)) then
      if (myrank == 0) call task_second_variation()
    end if

    if (allocated(evalfv)) deallocate(evalfv)
    if (allocated(occfv)) deallocate(occfv)
    call delete_selfenergy

    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)

    return
end subroutine
!EOC
