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
    call calcvxcnn
    call timesec(t1)

    ! clean not used anymore global exciting variables
    call clean_gndstate

    if (input%gw%taskname /= 'g0w0-x') then
      if (.not.input%gw%rpmat) then
        !========================================================
        ! calculate momentum matrix elements and store to a file
        !========================================================
        call calcpmatgw
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

    if (myrank==0) call boxmsg(fgw,'=','GW cycle')

    ! each process does a subset
    do iq = iqstart, iqend

      if (myrank==0) then
        write(fgw,*) '(task_gw): q-point cycle, iq = ', iq
        call flushifc(fgw)
      end if

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

      if (input%gw%taskname /= 'g0w0-x') then
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
        call init_dielectric_function(mbsiz, iomstart, iomend, Gamma)
        select case (trim(input%gw%scrcoul%scrtype))
          case('ppm','PPM')
            call calcepsilon_ppm(iq, iomstart, iomend)
          case default
            call calcepsilon(iq, iomstart, iomend)
            !==========================================
            ! Calculate the screened Coulomb potential
            !==========================================
            call calcinveps(iomstart, iomend)
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
    if ((nproc_row>1) .and. (myrank_col==0)) then
      call mpi_sum_array(0, selfex, nbandsgw, kset%nkpt, mycomm_row)
      if (input%gw%taskname /= 'g0w0-x') then
        ! G0W0
        call mpi_sum_array(0, selfec, nbandsgw, freq_selfc%nomeg, kset%nkpt, mycomm_row)
        if (input%gw%taskname == 'cohsex') then
          call mpi_sum_array(0, sigsx, nbandsgw, kset%nkpt, mycomm_row)
          call mpi_sum_array(0, sigch, nbandsgw, kset%nkpt, mycomm_row)
        end if
      end if ! selfec
    endif
#endif

    !===============================================================================
    ! output block
    !===============================================================================

    if (myrank == 0) then

      if ((input%gw%taskname /= 'g0w0-x') .and. (input%gw%selfenergy%method == "ac")) then
        ! Analytical continuation of the correlation self-energy from the complex to the real frequency axis
        if (input%gw%printSelfC) call plot_selfc_iw()
        call calcselfc_ac()
      end if

!$OMP critical

      !===============================
      ! Write self-energies to files
      !===============================
      call write_selfenergy(ibgw, nbgw, kset%nkpt, freq_selfc%nomeg)

      !=======================================
      ! Calculate the quasiparticle energies
      !=======================================

      ! KS band structure
      evalks(ibgw:nbgw,:) = evalfv(ibgw:nbgw,:)
      ! call bandstructure_analysis('KS', &
      !     ibgw, nbgw, kset%nkpt, evalks(ibgw:nbgw,:), efermi)

      ! solve QP equation
      call calcevalqp()
      call write_qp_energies('EVALQP.DAT')
      call putevalqp('EVALQP.OUT', kset, ibgw, nbgw, evalks, eferks, evalqp, eferqp)

      if (.not.isspinorb()) then

        if (input%gw%taskname /= 'g0w0-x') then
          if (input%gw%printSelfC)            call plot_selfc()
          if (input%gw%printSpectralFunction) call plot_spectral_function()
        end if

        ! G0W0 QP band structure
        select case (input%gw%taskname)

          case('g0w0-x')
            call bandstructure_analysis('G0W0-X band structure', &
                ibgw, nbgw, kset%nkpt, evalqp(ibgw:nbgw,:), eferqp)

          case('cohsex')
            call bandstructure_analysis('COHSEX band structure', &
                ibgw, nbgw, kset%nkpt, evalqp(ibgw:nbgw,:), eferqp)

          case('g0w0')
            call bandstructure_analysis('G0W0 band structure', &
                ibgw, nbgw, kset%nkpt, evalqp(ibgw:nbgw,:), eferqp)

        end select

      end if

!$OMP end critical

    end if ! myrank
    call barrier() ! synchronize all threads

    !-----------------------------------------
    ! Second-variational treatment of SO
    !-----------------------------------------
    if (isspinorb()) then
      call init0()
      call readstate()
      if (myrank==0) call task_second_variation()
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
