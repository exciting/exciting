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
    use calculate_correlation_self_energy, only: calcselfc, sigmac_indexes
    use calculate_dielectric_function, only: calcepsilon, epsilon_indexes
    use calculate_exchange_self_energy, only: calcselfx
    use constants, only: zzero
    use invert_dielectric_function, only: calcinveps
    use mod_bands, only: bandstructure_analysis, delete_bands, evalfv, nstdf, nstse, numin, occfv
    use mod_core_states, only: n_core_states => ncg
    use mod_coulomb_potential, only: barc, delete_coulomb_potential, calculate_singularities_coeff
    use mod_dielectric_function, only: eps00, epsh, epsw1, epsw2, epsilon, init_dielectric_function, &
      delete_dielectric_function
    use mod_frequency, only: delete_freqgrid
    use mod_gw_degeneracies, only: ibgw_including_degeneracy, nbgw_including_degeneracy
    use mod_gaunt_coefficients, only: delete_gaunt_coefficients
    use mod_kpointset, only: delete_Gk_vectors, delete_k_vectors, delete_kq_vectors, delete_G_vectors
    use mod_mpi_gw, only: iomstart, iomend, iqstart, iqend, mpi_sum_array, indexes_parallelization
    use mod_misc_gw, only: Gamma, gammapoint
    use mod_product_basis, only: mpwipw, locmatsiz, mbsiz, matsiz, delete_product_basis
    use mod_selfenergy, only: evalks, evalqp, eferks, eferqp, znorm, singc1, singc2, &
      selfec, selfex, freq_selfc, sigc, sigsx, sigch, plot_selfc, plot_selfc_iw, &
      init_selfenergy, write_selfenergy_binary, delete_selfenergy
    use mod_vxc, only: calcvxcnn, write_vxcnn, vxcnn
    use modinput, only: input, isspinorb
    use modmpi, only: barrier, distribute_loop, mpiglobal, rank
    use modgw, only: fgw, kset, kqset, Gqset, Gkset, Gkqset, Gset, Gqbarc, ibgw, nbgw, nbandsgw, &
      ciw, kiw, unw, kcw, freq, time_dfinv
    use modxs, only: symt2
    use quasiparticle_energies, only: write_qp_energies_text_format
    use mod_APW_LO, only: lorbl, nlorb, apword
    use mod_atoms, only: idxas
    use mod_eigensystem, only: idxlo
    use mod_muffin_tin, only: idxlm
    use precision, only: dp, i32
#include "offload.fpp"
    
!!LOCAL VARIABLES:
    implicit none

    integer(i32) :: iq, ik, mdim

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
    if (input%gw%coreflag=='all') then
      mdim = nstse+n_core_states
    else
      mdim = nstse
    end if

    !=================================================
    ! Calculate the diagonal matrix elements of the
    ! DFT exchange-correlation potential
    !=================================================
    ! it is better to do it here to deallocate cfunir and vxcir arrays
    call calcvxcnn( ibgw_including_degeneracy, nbgw_including_degeneracy, [(ik, ik=1,kset%nkpt)], kset%vkl(:, 1:kset%nkpt), mpiglobal )
    if (rank==0) then
      call write_vxcnn( 'binary', ibgw, nbgw )
      call write_vxcnn( 'text', ibgw, nbgw )
    end if

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
    singc1 = 0.d0
    singc2 = 0.d0
  
    !---------------------------------------
    ! treatment of singularities at G+q->0
    !---------------------------------------
    call calculate_singularities_coeff( input%gw%barecoul%cutofftype, &
      input%gw%selfenergy%singularity, kqset%nkpt, singc2 )

    ! initialize self-energy arrays
    call init_selfenergy(ibgw, nbgw, kset%nkpt)

    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================
    call distribute_loop( mpiglobal, kqset%nkpt, iqstart, iqend )
    iomstart = 1
    iomend = freq%nomeg

    if (rank==0) call boxmsg(fgw,'=','GW cycle')

    ! each process does a subset
    do iq = iqstart, iqend

      if (rank==0) then
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
      call calcselfx( iq, 1, kset%nkpt )

      if (input%gw%taskname /= 'g0w0-x') then
        !========================================
        ! Set v-diagonal MB and reduce its size
        !========================================
        call setbarcev(input%gw%barecoul%barcevtol, Gamma)
        call delete_coulomb_potential
        !===================================
        ! Calculate the dielectric function
        !===================================
        call init_dielectric_function(mbsiz, iomstart, iomend, Gamma)
        select case (trim(input%gw%scrcoul%scrtype))
          case('ppm','PPM')
            call calcepsilon_ppm(iq, iomstart, iomend)
          case default
            call calcepsilon(iq, epsilon_indexes( &
                            indexes_parallelization( 1, kqset%nkpt, 1, kqset%nkpt ), &
                            indexes_parallelization( numin, nstdf, numin, nstdf ), &
                            indexes_parallelization( iomstart, iomend, iomstart, iomend ) ) &
                            )
            !==========================================
            ! Calculate the screened Coulomb potential
            !==========================================
            if( gamma ) then 
              call calcinveps(iomstart, iomend, gamma, input%gw%scrcoul, freq%fconv, symt2,&
                              &epsilon, epsw1, epsw2, epsh, eps00, time_dfinv)
            else
              call calcinveps( iomstart, iomend, gamma, freqtype=freq%fconv, epsilon=epsilon, time_dfinv=time_dfinv )
            end if 
        end select
        !========================================
        ! Calculate the q-dependent self-energy
        !========================================
        call calcselfc( iq, sigmac_indexes( indexes_parallelization( 1, kset%nkpt, 1, kset%nkpt ), &
                                            indexes_parallelization( 1, mdim, 1, mdim ) &
                                          ) &
                      )
        call delete_dielectric_function(Gamma)
        if (allocated(kcw)) deallocate(kcw)
        if (allocated(unw)) deallocate(unw)
      end if

      ! clean unused data
      if (allocated(mpwipw)) then
        OMP_OFFLOAD target exit data map(delete: mpwipw)
        deallocate(mpwipw)
      end if
      if (allocated(barc)) then
          OMP_OFFLOAD target exit data map(delete: barc)
          deallocate(barc)
      end if
      !call omp_set_num_threads(nthreads)
    end do ! iq

    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)

    call mpi_sum_array( selfex, mpiglobal, .false. )
    if (input%gw%taskname /= 'g0w0-x') then
      ! G0W0
      call mpi_sum_array( selfec, mpiglobal, .false. )
      if (input%gw%taskname == 'cohsex') then
        call mpi_sum_array( sigsx, mpiglobal, .false. )
        call mpi_sum_array( sigch, mpiglobal, .false. )
      end if
    end if ! selfec

    !===============================================================================
    ! output block
    !===============================================================================

    if (rank == 0) then
      if ((input%gw%taskname /= 'g0w0-x') .and. (input%gw%selfenergy%method == "ac")) then
        ! Analytical continuation of the correlation self-energy from the complex to the real frequency axis
        if (input%gw%printSelfC) call plot_selfc_iw()
        call calcselfc_ac()
      end if

!$OMP critical

      !===============================
      ! Write self-energies to files
      !===============================
      call write_selfenergy_binary(ibgw, nbgw, kset%nkpt, freq_selfc%nomeg)

      !=======================================
      ! Calculate the quasiparticle energies
      !=======================================

      ! KS band structure
      evalks(ibgw:nbgw,:) = evalfv(ibgw:nbgw,:)

      ! solve QP equation
      call calcevalqp()
      ! Write QP energies into an output file
      select case(input%gw%taskname)
        case('g0w0')
          call write_qp_energies_text_format( [(ik, ik=1,kset%nkpt)], kset%vkl(:, 1:kset%nkpt), kset%wkpt, &
            ibgw, evalks, evalqp, real( vxcnn%diag_elements(ibgw:, :), dp ), selfex, sigc, znorm )
        
        case('g0w0-x')
          if( allocated(sigc) ) deallocate( sigc )
          allocate( sigc(ibgw:nbgw, kset%nkpt), source=zzero )
          if( allocated(znorm) ) deallocate( znorm )
          allocate( znorm(ibgw:nbgw, kset%nkpt), source=0._dp )
          call write_qp_energies_text_format( [(ik, ik=1,kset%nkpt)], kset%vkl(:, 1:kset%nkpt), kset%wkpt, &
            ibgw, evalks, evalqp, real( vxcnn%diag_elements(ibgw:, :), dp ), selfex, sigc, znorm )

        case('cohsex')
          if( allocated(znorm) ) deallocate( znorm )
          allocate( znorm(ibgw:nbgw, kset%nkpt), source=0._dp )
          call write_qp_energies_text_format( [(ik, ik=1,kset%nkpt)], kset%vkl(:, 1:kset%nkpt), kset%wkpt, &
            ibgw, evalks, evalqp, real( vxcnn%diag_elements(ibgw:, :), dp ), sigsx, sigch, znorm )
      
      end select
      call putevalqp('EVALQP.OUT', kset, ibgw, nbgw, evalks, eferks, evalqp, eferqp)

      if (.not.isspinorb()) then

        if (input%gw%taskname /= 'g0w0-x') then
          if (input%gw%printSelfC)            call plot_selfc(freq_selfc%freqs, [(ik, ik=1,kset%nkpt)], selfec, first_band=1)
          if (input%gw%printSpectralFunction) call plot_spectral_function()
        end if

        ! G0W0 QP band structure
        select case (input%gw%taskname)

          case('g0w0-x')
            call bandstructure_analysis('G0W0-X band structure', &
                ibgw, evalqp(ibgw:nbgw,:), eferqp, .true.)

          case('cohsex')
            call bandstructure_analysis('COHSEX band structure', &
                ibgw, evalqp(ibgw:nbgw,:), eferqp, .true.)

          case('g0w0')
            call bandstructure_analysis('G0W0 band structure', &
                ibgw, evalqp(ibgw:nbgw,:), eferqp, .true.)

        end select

      end if

!$OMP end critical

    end if ! rank
    call barrier() ! synchronize all threads

    !-----------------------------------------
    ! Second-variational treatment of SO
    !-----------------------------------------
    if (isspinorb()) then
      call init0()
      call readstate()
      if (rank==0) call task_second_variation()
    end if

    if (allocated(evalfv)) deallocate(evalfv)
    if (allocated(occfv)) deallocate(occfv)
    call delete_selfenergy

    OMP_OFFLOAD target exit data map(delete: kset, Gset, Gkset, Gkqset, Gqset, Gqbarc, kqset) 

    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_Gk_vectors(Gkqset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    call delete_bands()
    call delete_gaunt_coefficients()
    call delete_product_basis()

    OMP_OFFLOAD target exit data map(delete: idxas, idxlo, idxlm, lorbl, apword, nlorb)

end subroutine
!EOC
