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
    
    !> Cubic GW (Manoar Hossain)
    use precision, only: wp
    use mod_frequency, only: minimax_grid_type 
    use mod_dielectric_function, only: pola_q_mtmt, pola_q_mti, pola_q_ii
    use mod_polarizability_q_mtmt, only: polarizability_q_mtmt, polarizability_tau_to_omega_mtmt
    ! use mod_polarizability_R_mti, only: polarizability_R_mti
    use mod_polarizability_q_mti, only: polarizability_q_mti, polarizability_tau_to_omega_mti
    ! use mod_green_R_ii, only: green_R_ii
    ! use mod_polarizability_R_ii, only: polarizability_R_ii
    use mod_polarizability_q_ii, only: polarizability_q_ii, polarizability_tau_to_omega_ii
    ! TODO: Dress with IFDEF 
    use gx_minimax, only: gx_minimax_grid, gx_get_error_message 
    !! TEST and DELETE 
    use modmain, only : ngkmax   !! max value of G-vectors in all k-points
  

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

    ! minimax grid
    !> Minimum and maximum transition energies
    real(dp) :: e_transition_min, e_transition_max
    type(minimax_grid_type) :: minimax_grid
    real(dp) :: gx_max_errors(3)
    real(dp) :: gx_cosft_duality_error
    integer  :: gx_ierr
    integer :: it, jt
    character(len=500) :: gx_msg

    ! Local variables for cubic GW (Manoar)
    integer :: ir, nbigR
    integer :: nlm_tot
    integer :: ias1, ias2
    integer :: nlm1, nlm2 
    integer :: igq1, igq2
    real(wp) :: t_i, t_f
    complex(wp), allocatable :: pola_q_mtmt_omega(:,:,:,:,:,:)
    complex(wp), allocatable :: pola_q_mti_omega(:,:,:,:,:)
    complex(wp), allocatable :: pola_q_ii_omega(:,:,:,:)


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

    !-------------------------------------------------------------------------------------------------
    !Cubic GW ---starts
    select case (input%gw%taskname)
      case ('g0w0_cubic')

      ! Initialise the imaginary times, imaginary frequencies and the transformation weights
      ! Note, `freq%nomeg` is a placeholder name (gives grid points)
      call minimax_grid%init('imag', input%gw%freqgrid%nomeg)

      
      !(nomax and numin) are indices of VBM and CBM, respectively
      write(*, *) evalfv(numin,:) 
      write(*, *) 
      write(*, *) evalfv(nomax,:) 
      write(*, *) 
      write(*, *) evalfv(:,1) 
      
      e_transition_min = minval(evalfv(numin,:)) - maxval(evalfv(nomax,:))
      e_transition_max = maxval(evalfv(nstfv,:)) - minval(evalfv(1,:)) 
      write(*, *) 'Minimum and maximum transition energies, respectively', e_transition_min, e_transition_max
      write(*, *) 'N freq', minimax_grid%n_points 

      call gx_minimax_grid(minimax_grid%n_points, e_transition_min, e_transition_max, & 
                           minimax_grid%time, minimax_grid%time_weights, & 
                           minimax_grid%freq, minimax_grid%freq_weights, & 
                           minimax_grid%cos_time_to_freq_weights, & 
                           minimax_grid%cos_freq_to_time_weights, & 
                           minimax_grid%sin_time_to_freq_weights, & 
                           gx_max_errors, gx_cosft_duality_error, gx_ierr) 

      if (gx_ierr /= 0) then 
        call gx_get_error_message(gx_msg) 
        write(*, *) trim(gx_msg) 
      endif 
      
      write(1, *) 'Img time grid n-points', minimax_grid%n_points
      write(1, *) 'tau                                omega'
      do it = 1 , minimax_grid%n_points
        write(1, *) minimax_grid%time(it), minimax_grid%freq(it)
      enddo
      do it = 1 , minimax_grid%n_points
        do jt = 1 , minimax_grid%n_points
          write(7,*) minimax_grid%cos_time_to_freq_weights(jt, it)
        enddo 
      enddo 
 
      ! Initialise variables, polarizability
      nlm_tot = locmatsiz/natmtot
      nbigR = kqset%nkpt
      !-------------------------------------------------------------------------------------------------
      !allocate(polarisability_mtmt(locmatsiz,locmatsiz,kqset%nkpt,minimax_grid%n_points))
      !print*, 'Polarizability: shape', shape(polarisability_mtmt), 'size =', size(polarisability_mtmt)
      !print*, 'freq grid', freq%nomeg, input%gw%freqgrid%nomeg
      print*, '-------------------------------------------------------------------------------------------------'
      print*, 'testPrint variables'
      print*, 'Total atoms', natmtot
      print*, 'Total k-points (irreduced) and Bravais latt.', nkpt, kset%nkptnr, kqset%nkpt, nbigR
      print*, 'Total k-points (reduced)', kset%nkpt, size(evalfv,2)
      print*, 'Fermi energy/chemical potential', efermi
      print*, 'VBM and CBM positions', nomax, numin
      print*, 'combined "lm" dimension ->lmmaxapw', lmmaxapw
      print*, 'Number of states used to calculate the dielectric function', nstdf
      print*, 'Number of states used to calculate the self-energy', nstse
      print*, 'Size of the local part of the mixed basis including LM combinations locmatsiz', locmatsiz
      print*, 'locmatsiz/natmtot, nlm_tot', locmatsiz/natmtot, nlm_tot
      !print*, 'integer k', ikv
      !print*, 'mbindex test', mbindex(:,4)!, mbindex(1,5)
      print*, '-------------------------------------------------------------------------------------------------'
      ! stop

      ! lmmaxapw -> combined "lm" dimension
      ! natmtot --> total number of atoms in the system
      ! kqset%nkpt -> Total(irreduced) k-points. 
      ! NOTE: later it shoule be replaced by number of bravias vector R
      allocate(pola_q_mtmt_omega(nlm_tot, natmtot, nlm_tot, natmtot, kqset%nkpt, minimax_grid%n_points))
      allocate(pola_q_mti_omega(nlm_tot, natmtot, ngkmax, kqset%nkpt, minimax_grid%n_points))
      allocate(pola_q_ii_omega(ngkmax, ngkmax, kqset%nkpt, minimax_grid%n_points))

      pola_q_mtmt_omega = cmplx(0.0_wp, 0.0_wp, kind=wp)
      pola_q_mti_omega = cmplx(0.0_wp, 0.0_wp, kind=wp)
      pola_q_ii_omega = cmplx(0.0_wp, 0.0_wp, kind=wp)

      call CPU_TIME(t_i)
      
      do it = 1, minimax_grid%n_points
        print*, '-------------------------------------------------------------------------------'
        print*, '=========================== *** tau loop ***', it, '==========================='
        print*, '-------------------------------------------------------------------------------'
        
        call polarizability_q_mtmt(minimax_grid%time(it), pola_q_mtmt)
        call polarizability_tau_to_omega_mtmt(it, &
             minimax_grid%cos_time_to_freq_weights, &
             pola_q_mtmt, pola_q_mtmt_omega) 


        call polarizability_q_mti(minimax_grid%time(it), pola_q_mti)
        call polarizability_tau_to_omega_mti(it, &
             minimax_grid%cos_time_to_freq_weights, &
             pola_q_mti, pola_q_mti_omega) 


        ! call polarizability_q_ii(minimax_grid%time(it), pola_q_ii)
        ! call polarizability_tau_to_omega_ii(it, &
        !      minimax_grid%cos_time_to_freq_weights, &
        !      pola_q_ii, pola_q_ii_omega) 

      enddo   !! it loop ends ($\tau$)

      !! Polarizability in omega MT-MT TEST ------------------
      do it = 1, minimax_grid%n_points
        do iq = 1, kqset%nkpt
          do ias2 = 1, natmtot
            do nlm2 = 1, nlm_tot
              do ias1 = 1, natmtot
                do nlm1 = 1, nlm_tot
                  ! write(90+it, *) pola_q_mtmt_omega(nlm1, ias1, nlm2, ias2, iq, it)
                  write(1000000+1000*iq+it, *) pola_q_mtmt_omega(nlm1, ias1, nlm2, ias2, iq, it)
                enddo 
              enddo 
            enddo 
          enddo 
        enddo 
      enddo 
      !! Polarizability in omega MT-MT TEST ------------------

      !! Polarizability in omega MT-I TEST ------------------
      do it = 1, minimax_grid%n_points
        do iq = 1, kqset%nkpt
          do igq2 = 1, ngkmax
            do ias1 = 1, natmtot
              do nlm1 = 1, nlm_tot
                write(2000000+1000*iq+it, *) pola_q_mti_omega(nlm1, ias1, igq2, iq, it)
                ! write(20000+100*iq+70+it, *) pola_q_mtmt_omega(nlm1, ias1, nlm2, ias2, iq, it)
              enddo 
            enddo 
          enddo 
        enddo 
      enddo 
      !! Polarizability in omega MT-I TEST ------------------

      !! Polarizability in omega I-I TEST ------------------
      do it = 1, minimax_grid%n_points
        do iq = 1, kqset%nkpt
          do igq1 = 1, ngkmax  !! need to change with actual dimension
            do igq2 = 1, ngkmax  !! need to change with actual dimension
              write(3000000+1000*iq+it, *) pola_q_ii_omega(igq2, igq1, iq, it)
            enddo 
          enddo 
        enddo 
      enddo 
      !! Polarizability in omega I-I TEST ------------------
      
      call CPU_TIME(t_f)
      print*, 'Total time taken (cubic GW):', t_f - t_i
      
      ! stop
     !Cubic GW ---ends
     !-------------------------------------------------------------------------------------------------
      

      
    case default ! quartic algorithm follows

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


    end select ! Cubic vs quartic GW


    return
end subroutine
!EOC
