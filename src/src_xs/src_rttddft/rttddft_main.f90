! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created Apr 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Cleaned May 2023 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> This module is the kernel of a RT-TDDFT calculation.
!> It contains the subroutine `coordinate_rttddft_calculation`, which manages
!> a RT-TDDFT calculation. Also here are implemented the following subroutines:
!> uppot, uprho.
module rttddft_main
  use asserts, only: assert
  use errors_warnings, only: terminate_if_false
  use m_getunit, only: getunit
  use MD, only: force, MD_input_keys
  use MD_io, only: MD_out
  use mod_atoms, only: natmtot, natoms, nspecies, atposc, idxas
  use mod_charge_and_moment, only: chgval
  use mod_kpoint, only: nkpt
  use mod_lattice, only: omega
  use mod_misc, only: filext
  use mod_mpi_env, only: mpiinfo
  use modinput, only: input, input_type
  use modmpi, only: rank, procs, mpi_env_k, mpiglobal, distribute_loop, barrier
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_CurrentDensity, only: UpdateCurrentDensity
  use rttddft_Energy, only: TotalEnergy, obtain_energy_rttddft
  use rttddft_GlobalVariables
  use rttddft_HamiltonianOverlap, only: UpdateHam
  use rttddft_init, only: initialize_rttddft
  use rttddft_io, only: open_files_jpa, close_files_jpa, write_jpa, &
    open_file_timing, close_file_timing, write_timing, &
    open_file_nexc, close_file_nexc, write_nexc, &
    open_file_etot, close_file_etot, write_total_energy, &
    open_file_info, close_file_info, write_file_info, write_file_info_header, &
    write_wavefunction
  use rttddft_MD, only: force_rttdft, move_ions, &
    MD_allocate_global_arrays => allocate_global_arrays, &
    MD_deallocate_global_arrays => deallocate_global_arrays, &
    MD_evaluate_charge_val => evaluate_charge_val
  use rttddft_NumberExcitations, only: Obtain_number_excitations
  use rttddft_screenshot, only: screenshot
  use rttddft_VectorPotential, only: Calculate_Vector_Potential, Evolve_A_ind => Solve_ODE_Vector_Potential
  use rttddft_Wavefunction, only: UpdateWavefunction, Update_basis_derivative
  
  implicit none

  private
  
  public  :: coordinate_rttddft_calculation

contains

  !> This subroutine manages a RT-TDDFT calculation.
  !> <ol>
  !> <li> Run a single-shot groundstate calculation using the already converged
  !> density and potential. </li>
  !> <li> Obtain the KS wavefunctions, the density and the hamiltonian at 
  !> \( t=0 \). </li>
  !> <li> Evolve the wavefunctions, the density and the hamiltonian using the 
  !> desired time step. </li>
  !> </ol>
  subroutine coordinate_rttddft_calculation()
    ! counter for the number of iterations of real-time
    integer                 :: it
    ! prints data every nprint steps of the counter "it"
    integer                 :: nprint
    ! indexes of the first and the last k-points
    integer                 :: first_kpt, last_kpt
    integer                 :: niter_screenshot

    integer                 :: iprint
    integer                 :: is, ia, ias, timeStepMultiplier
    logical                 :: predCorrReachedMaxSteps

    character(len=100)      :: string
    character(len=100)      :: fieldType
    character(len=100)      :: vectorPotentialSolver

    real(dp), allocatable   :: atom_positions(:, :) ! in cartesian coordinates x, y, z
    real(dp), allocatable   :: atom_velocities(:, :) ! in cartesian coordinates x, y, z
    type(force)             :: forces
    type(MD_input_keys)     :: molecular_dynamics

    real(dp),allocatable    :: nex(:), ngs(:), nt(:)
    real(dp)                :: aindsave(3),pvecsave(3)
    real(dp)                :: jindsave(3),aextsave(3),atotsave(3)
    real(dp)                :: timei, timef, timesave, timeiter, timeaux
    type(MD_out)            :: MD_outputs

    ! Variables to store data and print
    real(dp),allocatable    :: timestore(:),aindstore(:,:),atotstore(:,:)
    real(dp),allocatable    :: jindstore(:,:), pvecstore(:,:)
    real(dp),allocatable    :: atposcstore(:,:,:), velstore(:,:,:)
    type(force),allocatable :: forces_store(:)
    logical,allocatable     :: printforces(:)
    logical                 :: take_screenshots
    logical, allocatable    :: screenshot_was_taken(:)

    type(TotalEnergy), allocatable  :: etotstore(:)
    type(Timing_RTTDDFT_and_MD)     :: timing
    type(Timing_RTTDDFT_and_MD), allocatable :: timing_store(:)


    call timesec(timesave)

    ! Interface with input parameters
    nprint = input%xs%realTimeTDDFT%printAfterIterations
    fieldType = input%xs%realTimeTDDFT%laser%fieldType
    vectorPotentialSolver = input%xs%realTimeTDDFT%vectorPotentialSolver
    take_screenshots = associated( input%xs%realTimeTDDFT%screenshots )
    if( take_screenshots ) niter_screenshot = input%xs%realTimeTDDFT%screenshots%niter

    call molecular_dynamics%parse_input()
    ! we only perform MD in RT-TDDFT if the type is Ehrenfest
    if( molecular_dynamics%on ) molecular_dynamics%on = ( trim(molecular_dynamics%MD_type) == 'Ehrenfest' )

    call sanity_checks( input, nprint, mpiglobal )
    
    ! Outputs general info to RTTDDFT_INFO.OUT and
    ! opens TIMING_RTTDDFT.OUT (if this is the case)
    if(input%xs%realTimeTDDFT%printTimingGeneral) call timesec (timei)
    if( rank == 0 ) then
      call open_file_info
      call write_file_info_header
    end if

    call initialize_rttddft( molecular_dynamics )
    
    if( molecular_dynamics%on ) call init_MD( tstep, timeStepMultiplier, molecular_dynamics, &
      MD_outputs, atom_positions, atom_velocities, forces )

    ! Allocate variables to be stored and printed only after nprint steps
    allocate(timestore(nprint), aindstore(3,nprint), atotstore(3,nprint))
    allocate(jindstore(3,nprint), pvecstore(3,nprint))
    if( printTimesGeneral ) then
      allocate( timing_store(nprint) )
      allocate( screenshot_was_taken(nprint), source=.False. )
    end if
    if( calculateTotalEnergy ) allocate(etotstore(nprint))
    if( calculateNexc ) allocate(nex(nprint),ngs(nprint),nt(nprint))
    if( molecular_dynamics%on ) then
      allocate( printforces(nprint), atposcstore(3,natmtot,nprint), velstore(3,natmtot,nprint))
      if ( molecular_dynamics%print_all_force_components ) then
        allocate( forces_store(nprint) )
        do is = 1, nprint
          call forces_store(is)%allocate_arrays( natmtot )
        end do
      end if
    end if

    if( rank == 0 ) then
      call open_files_jpa
      call write_jpa( time, aind, atot, label='avec' )
      call write_jpa( time, pvec, label='pvec' )
      call write_jpa( time, jind, label='jind' )
    end if

    ! Initialize integers that contain the first and last k-point
    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

    ! Total energy
    if ( calculateTotalEnergy ) then
      call potcoul
      call potxc
      call obtain_energy_rttddft( first_kpt, last_kpt, ham_time, evecfv_gnd, etotstore(1) )
      ! Trick: we need an array to call the subroutine print_total_energy
      timestore(1) = time
      if( rank == 0 ) then
        call open_file_etot
        call write_total_energy( .True., 1, timestore(1), etotstore(1) )
      end if
    end if

    ! Number of excitations
    if (calculateNexc) then
      call Obtain_number_excitations( first_kpt, last_kpt, evecfv_gnd, &
        & evecfv_time, overlap, nex(1), ngs(1), nt(1) )
      ! Trick: we need an array to call the subroutine print_nexc
      timestore(1) = time
      if( rank == 0 ) then
        call open_file_nexc
        call write_nexc( .True., 1, timestore(1), nex(1), ngs(1), nt(1) )
      end if
    end if

    if ( take_screenshots ) call screenshot( 0, first_kpt, last_kpt, overlap, evecfv_gnd, &
        & evecfv_time, ham_time )

    if( printTimesGeneral ) then
      call timesec( timef )
      if( rank == 0 ) then 
        call open_file_timing
        call write_timing( timef-timesave ) !write time for initialization
      end if
    end if

    iprint = 1
    timeiter = timef
    ! This is the most important loop (performed for each time step \(\Delta t\)
    do it = 1, nsteps
      ! Variable to store the timing of each iteration
      timei = timeiter

      ! The "real time" t of our evolution
      time = time + tstep

      ! WAVEFUNCTION
      if ( predictorCorrector ) evecfv_save(:,:,:) = evecfv_time(:,:,:)
      if ( molecular_dynamics%on .and. molecular_dynamics%basis_derivative ) then
        call UpdateWavefunction( .False., atom_velocities )
      else 
        call UpdateWavefunction( .False. )
      end if
      if ( printTimesGeneral ) call timesecRTTDDFT(timei,timef,timing%t_RTTDDFT%t_wvf)

      ! Update the paramagnetic component of the induced current density
      call UpdateCurrentDensity( first_kpt, last_kpt, evecfv_time(:,:,:),jparanext(:) )
      if ( input%xs%realTimeTDDFT%subtractJ0 ) jparanext(:) = jparanext(:)-jparaspurious(:)
      if ( printTimesGeneral ) call timesecRTTDDFT(timei,timef,timing%t_RTTDDFT%t_curr)

      ! DENSITY
      if( printTimesGeneral ) then
        call uprho(it,timei,timef,timing%t_RTTDDFT)
      else 
        call uprho(it)
      end if

      ! KS-POTENTIAL
      if( printTimesGeneral ) then
        call uppot(timei,timef,timing%t_RTTDDFT)
      else 
        call uppot
      end if

      ! VECTOR POTENTIAL
      if(printTimesGeneral) call timesec(timei)
      ! Check if we need to save aind, pvec, atot and aext
      if( fieldType == 'external' .and. predictorCorrector .and. (vectorPotentialSolver /= 'euler') ) then
        aindsave(:) = aind(:)
        pvecsave(:) = pvec(:)
        aextsave(:) = aext(:)
      end if
      atotsave = atot
      if( fieldType == 'total' ) then
        call update_vector_potential( time, atot )
        if( molecular_dynamics%on ) then
          call Calculate_Vector_Potential( time+tstep, aindsave ) ! trick: aindsave is an auxiliary variable
          efield = obtain_electric_field( 2*tstep, aindsave, atotsave )
        end if
      else 
        call update_vector_potential( time, atot, aind, aext )
        if( molecular_dynamics%on ) efield = obtain_electric_field( tstep, atot, atotsave )
      end if
      if(printTimesGeneral) call timesecRTTDDFT(timei,timef,timing_store(iprint)%t_RTTDDFT%t_obtaina)

      ! INDUCED CURRENT
      ! Update the diamagnetic component of the induced current density
      jdia(:) = - atot(:) * chgval / c / omega
      ! Update the paramagnetic component of the induced current density
      jparaold(:) = jpara(:)
      jpara(:) = jparanext(:)
      ! Update the total induced current
      if ( predictorCorrector .and. (vectorPotentialSolver /= 'euler') ) then
        jindsave(:) = jind(:)
      end if
      jind(:) = jpara(:) + jdia(:)

      ! HAMILTONIAN
      call UpdateHam( predcorr=.False., calculateOverlap=.False., &
        timeGen=printTimesGeneral, timeDetail=printTimesDetailed, &
        timeini=timei, timefinal=timef, thmlint=timing%t_RTTDDFT%t_hmlint, &
        tham=timing%t_RTTDDFT%t_ham, &
        update_mathcalH=.False., update_mathcalB=.False., update_pmat=.False. )
      if( printTimesGeneral ) then
        timing%t_RTTDDFT%t_upham = timef - timei
        timei = timef
      end if
      ! Remark: it makes no sense to employ the predictor-corrector method with SE or EH!
      if ( predictorCorrector .and. (method /= 'SE') .and. (method /= 'EH') ) then
        call loopPredictorCorrector( it, maxstepsPredictorCorrector, &
          (fieldType == 'external') .and. (vectorPotentialSolver /= 'euler'), &
          first_kpt, last_kpt, aindsave, atotsave, aextsave, pvecsave, jindsave, &
          jparaold, predCorrReachedMaxSteps )
        
        if ( predCorrReachedMaxSteps .and. rank == 0 ) &
          write(*,*) 'Problems with convergence (PredCorr), time: ', time
        if ( molecular_dynamics%on .and. (fieldType == 'external')) &
          efield(:) = (-1d0/c/tstep)*(atot(:)-atotsave(:))
        if (printTimesGeneral) call timesecRTTDDFT(timei,timef,timing%t_RTTDDFT%t_predcorr)
      end if !predictor-corrector

      ! Obtain the total energy, if requested
      if( calculateTotalEnergy ) then
        call obtain_energy_rttddft( first_kpt, last_kpt, ham_time, evecfv_time, etotstore(iprint) )
        if ( printTimesDetailed ) call timesecRTTDDFT( timei, timef, timing%t_RTTDDFT%t_toten )
      end if

      ! Obtain the number of excited electrons, if requested
      if( calculateNexc ) then
        call Obtain_number_excitations( first_kpt, last_kpt, evecfv_gnd, &
          & evecfv_time, overlap, nex(iprint), ngs(iprint), nt(iprint))
        if( printTimesDetailed ) call timesecRTTDDFT( timei, timef, timing%t_RTTDDFT%t_nexc )
      end if

      if ( molecular_dynamics%on ) then
        if ( mod( it, timeStepMultiplier ) == 0 ) then
          if ( printTimesGeneral ) then 
            call timesec( timei )
            timing%t_Ehrenfest%MD_was_carried_out = .True.
          end if
          call forces%save_total_force()
          call force_rttdft( forces, molecular_dynamics%core_corrections, molecular_dynamics%valence_corrections,&
            printTimesDetailed, timei, timef, timing%t_Ehrenfest%t_MD_1st, &
            timing%t_Ehrenfest%t_MD_2nd, timing%t_Ehrenfest%t_MD_sumforces )
          call move_ions( forces%total, forces%total_save, molecular_dynamics%time_step, &
            atom_velocities, printTimesGeneral, printTimesDetailed, timei, timef, &
            timing%t_Ehrenfest%t_MD_moveions, timing%t_Ehrenfest%t_MD_updateBasis )
          printforces(iprint) = .True.
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia,is)
              atposcstore(1:3, ias, iprint) = atposc(1:3, ia, is)
            end do
          end do
          velstore(:,:,iprint) = atom_velocities(:,:)
          forces_store(iprint) = forces
          ! Update Hamiltonian with the new basis
          if( molecular_dynamics%update_overlap .or. allocated(mathcalH) .or. &
            & allocated(mathcalB) .or. molecular_dynamics%update_pmat ) then
            if(printTimesGeneral) call timesec(timeaux)
            call UpdateHam( predcorr=.False., &
              & calculateOverlap=molecular_dynamics%update_overlap, &
              & update_mathcalH=allocated(mathcalH), &
              & update_mathcalB=allocated(mathcalB), &
              & update_pmat=molecular_dynamics%update_pmat, &
                timeGen=printTimesGeneral, timeDetail=printTimesDetailed, &
                timeini=timeaux, timefinal=timef, tgenpmatbasis=timing%t_Ehrenfest%t_MD_pmat,&
                thmlint=timesave, tham=timing%t_Ehrenfest%t_MD_hamoverl )
            if( printTimesDetailed ) timing%t_Ehrenfest%t_MD_hamoverl = &
              timing%t_Ehrenfest%t_MD_hamoverl + timesave
          end if
          if( printTimesGeneral ) then
            timing%t_Ehrenfest%t_MD_step = timef - timei
            timei = timef
          end if
        else ! if ( mod( it, timeStepMultiplier ) == 0 )
          printforces(iprint) = .False.
          if ( printTimesGeneral ) timing%t_Ehrenfest%MD_was_carried_out = .False.
        end if ! if ( mod( it, timeStepMultiplier ) == 0 )
      end if ! if ( molecular_dynamics%on ) then

      ! Check if a screenshot has been requested
      if ( take_screenshots ) then
        if ( mod( it, niter_screenshot ) == 0 ) then
          if( printTimesGeneral ) screenshot_was_taken(iprint) = .True.
          if( printTimesGeneral ) call timesec(timei)
          call screenshot( it, first_kpt, last_kpt, overlap, evecfv_gnd, &
            & evecfv_time, ham_time )
          if( printTimesGeneral ) call timesecRTTDDFT( timei, timef, timing%t_RTTDDFT%t_screenshot )
        else 
          if( printTimesGeneral ) screenshot_was_taken(iprint) = .False.
        end if
      end if

      ! Store relevant information from this iteration
      timestore(iprint) = time
      aindstore(:, iprint) = aind(:)
      atotstore(:, iprint) = atot(:)
      pvecstore(:, iprint) = pvec(:)
      jindstore(:, iprint) = jind(:)
      if( printTimesGeneral ) timing_store(iprint) = timing

      ! Print relevant information, every 'nprint' steps
      if ( iprint == nprint ) then
        ! Update the counter
        iprint = 1
        if( rank == 0 ) then
          call write_jpa( timestore, aindstore, atotstore, label='avec' )
          call write_jpa( timestore, pvecstore, label='pvec' )
          call write_jpa( timestore, jindstore, label='jind' )
          if ( calculateTotalEnergy ) call write_total_energy( .False., nprint, &
            timestore(:), etotstore(:) )
          if ( calculateNexc ) call write_nexc( .False., nprint, timestore(:), &
            nex(:), ngs(:), nt(:) )

          ! Print forces - if this has been requested
          if( molecular_dynamics%on ) then
            do iprint = 1, nprint
              if( printforces(iprint) ) call write_MD_outputs( timestore(iprint), &
                atposcstore(:, :, iprint), velstore(:,:,iprint), forces_store(iprint), &
                molecular_dynamics%print_all_force_components, MD_outputs )
            end do
          end if ! if( molecular_dynamics%on )
        end if

        ! Update the counter
        iprint = 1
        if( printTimesGeneral ) then
          call timesecRTTDDFT( timeiter, timef, timing_store(nprint)%t_iteration )
          call write_timing( it, nprint, timing_store, screenshot_was_taken, molecular_dynamics%on )
        end if
      else ! if ( iprint .eq. nprint ) then
        if(printTimesGeneral) call timesecRTTDDFT(timeiter,timef,timing_store(iprint)%t_iteration)
        iprint = iprint + 1
      end if ! if ( iprint == nprint ) 
      ! Make all the processes wait here: the master alone has been writing the files above
      call barrier( mpi_env_k )
    end do ! do it = 1, nsteps

    if ( rank == 0 ) then
      call close_files_jpa
      call write_file_info( 'Real-time TDDFT calculation finished' )
      call close_file_info
      if( calculateTotalEnergy ) call close_file_etot
      if( calculateNexc ) call close_file_nexc
      if( printTimesGeneral ) call close_file_timing
      if( molecular_dynamics%on ) call MD_outputs%close_files()
    end if

    ! write wavefunction, and potential and density with _RTTDDFT.OUT as suffix
    string = filext
    filext = '_RTTDDFT'//trim(filext)
    call write_wavefunction( first_kpt, evecfv_time )
    if ( rank == 0 ) call writestate
    filext = string

    call deallocate_global_arrays( molecular_dynamics%on )
    
  end subroutine coordinate_rttddft_calculation

  !> This is just an interface to call the subroutine `[[UpdateDensity]]`, which
  !> updates the charge density
  subroutine uprho( iteration_counter, timei, timef, timing )
    use rttddft_Density, only: UpdateDensity
    !> Tells how many time steps have already been executed
    integer, intent(in) :: iteration_counter
    !> time (in seconds) when the subroutine was called (used for making time differences)
    real(dp), intent(inout), optional :: timei
    !> time (in seconds) after executing this subroutine
    real(dp), intent(inout), optional :: timef
    type(TimingRTTDDFT), optional     :: timing

    if( printTimesGeneral ) call assert( present(timei) .and. present(timef) &
      .and. present(timing), 'Optional arguments must be present if printTimesGeneral is true')

    if ( printTimesGeneral ) then
      if ( printTimesDetailed ) then
        call UpdateDensity( iteration_counter, timei, timef, timing%t_dens_rho, &
          & timing%t_dens_symrf, timing%t_dens_rfmtctof, &
          & timing%t_dens_addrhocr, timing%t_dens_charge, &
          & timing%t_dens_rhonorm )
      else
        call UpdateDensity( iteration_counter, timei, timef )
      end if
    else
      call UpdateDensity( iteration_counter )
    endif
    if ( printTimesGeneral ) then
      timing%t_dens = timef-timei
      timei = timef
    endif
  end subroutine uprho


  !> This is just an interface to call the subroutines that updates the KS potential
  subroutine uppot( timei, timef, timing )
    real(dp), intent(inout), optional :: timei
    real(dp), intent(inout), optional :: timef
    type(TimingRTTDDFT), optional     :: timing

    real(dp)             :: timeisave

    if( printTimesGeneral ) then 
      call assert( present(timei) .and. present(timef) &
      .and. present(timing), 'Optional arguments must be present if printTimesGeneral is true')
      timeisave = timei
    end if 

    ! Compute the effective potential (with the updated density)
    call poteff
    if( printTimesDetailed ) call timesecRTTDDFT( timei, timef, timing%t_poteff )
    ! Fourier transform effective potential to G-space
    call genveffig
    if( printTimesDetailed ) call timesecRTTDDFT( timei, timef, timing%t_genveffig )
    call genmeffig
    if( printTimesGeneral ) then
      call timesecRTTDDFT( timeisave, timef, timing%t_uppot )
      if( printTimesDetailed ) timing%t_genmeffig = timef-timei
      timei = timef
    end if
  end subroutine uppot

  !> (private subroutine) Check if variables given in the input file make sense
  subroutine sanity_checks( inp, nprint, mpi_env  )
    !> type with the variables given in the input file
    type(input_type):: inp
    !> Variable that tells after how many iterations the output should be written
    integer(i32), intent(in) :: nprint
    !> Variable that encapsulates the MPI communicator that will abort if a sanity test fails
    type(mpiinfo), intent(inout):: mpi_env
    
    call terminate_if_false( mpi_env, nprint > 0, &
      & 'Invalid non-positive number as printAfterIterations' )

    call terminate_if_false( mpi_env, .not. inp%groundstate%solver%packedmatrixstorage, &
      & 'Error: RT-TDDFT does not work with matrices stored in a packed form.' )

    ! Consistency check: check if no spin polarized calculations are requested.
    call terminate_if_false( mpi_env, .not. inp%groundstate%tevecsv, &
      & 'Error: only spin unpolarised calculations are possible with RT-TDDFT now.' )

    ! Consistency check: laser has been defined?
    call terminate_if_false( mpi_env, associated(inp%xs%realTimeTDDFT%laser), &
      & 'Element <laser> in <realTimeTDDFT> not found')

    ! Consistency check: MD and predictor corrector?
    call terminate_if_false( mpi_env, .not. ( associated(inp%xs%realTimeTDDFT%predictorCorrector) &
      .and. associated(inp%MD) ), &
      & 'It is currently not possible to use the predictor corrector method together with molecular dynamics')

  end subroutine

  subroutine update_vector_potential( t, A_tot, A_ind, A_ext )
    !> time \(t\)
    real(dp), intent(in)              :: t
    !> total vector potential: \(A_{tot} = A_{ind} + A_{ext}\)
    real(dp), intent(inout)           :: A_tot(3)
    !> induced vector potential
    real(dp), intent(inout), optional :: A_ind(3)
    !> external vector potential
    real(dp), intent(inout), optional :: A_ext(3)

    logical :: all_fields_present

    all_fields_present = present(A_ind)
    if( all_fields_present ) call assert( present(A_ext), 'If A_ind is present, then A_ext must also be' )

    !TODO(Ronaldo): Refactor `Evolve_A_ind` to avoid globals
    call Evolve_A_ind( .not. all_fields_present )
    if( all_fields_present ) then
      call Calculate_Vector_Potential( t, A_ext(:) )
      A_tot(:) = A_ind(:) + A_ext(:)
    else
      call Calculate_Vector_Potential( t, A_tot(:) )
    end if
  end subroutine

  !> Obtain the electric field as the time derivative of the vector potential
  pure function obtain_electric_field( dt, A_tot, A_tot_previous ) result( E_field )
    !> time step
    real(dp), intent(in)  :: dt
    !> vector potential at time `t`
    real(dp), intent(in)  :: A_tot(3)
    !> vector potential at time `t-dt`
    real(dp), intent(in)  :: A_tot_previous(3)
    real(dp) :: E_field(3)
    E_field = ( -1._dp / c / dt ) * ( A_tot - A_tot_previous )
  end function

  subroutine loopPredictorCorrector( it, maxSteps, evolveA, first_kpt, last_kpt, &
    aindsave, atotsave, aextsave, pvecsave, jindsave, jparasave, maxStepsReached )
    !> current iteration number in the RT-TDDFT loop
    integer(i32), intent(in)       :: it
    !> maximum steps in the predictor corrector scheme
    integer(i32), intent(in)       :: maxSteps
    !> If `.True.`, the vector potential is evolved in each step
    logical, intent(in)            :: evolveA
    !> index of the first `k-point` to be considered in the sum
    integer(i32),intent(in)        :: first_kpt
    !> index of the last `k-point` considered
    integer(i32),intent(in)        :: last_kpt
    !> Backup of `aind`
    real(dp), intent(in)           :: aindsave(3)
    !> Backup of `atot`
    real(dp), intent(in)           :: atotsave(3)
    !> Backup of `aext`
    real(dp), intent(in)           :: aextsave(3)
    !> Backup of `pvec`
    real(dp), intent(in)           :: pvecsave(3)
    !> Backup of `jind`
    real(dp), intent(in)           :: jindsave(3)
    !> Backup of `jpara`
    real(dp), intent(in)           :: jparasave(3)
    !> When `.True.`, it informs that the maximum steps have been reached
    logical, intent(out)           :: maxStepsReached

    integer(i32) :: i
    real(dp)     :: err

    do i = 1, maxSteps
      ! WAVEFUNCTION
      evecfv_time(:,:,:) = evecfv_save(:,:,:)
      call UpdateWavefunction( .True. )

      ! Update the paramagnetic component of the induced current density
      call UpdateCurrentDensity( first_kpt, last_kpt, evecfv_time(:,:,:),jparanext(:))
      if ( input%xs%realTimeTDDFT%subtractJ0 ) jparanext(:) = jparanext(:)-jparaspurious(:)

      ! DENSITY
      call uprho( it )

      ! KS-POTENTIAL
      call uppot

      ! VECTOR POTENTIAL
      ! Update the induced part of the vector potential
      if( evolveA ) then
        jpara(:) = jparasave(:) !attention: jparaold saves the value of jpara(t-deltat)
        aind(:) = aindsave(:)
        atot(:) = atotsave(:)
        aext(:) = aextsave(:)
        pvec(:) = pvecsave(:)
        jind(:) = jindsave(:)
        call update_vector_potential( time, atot, aind, aext )
        call Evolve_A_ind( .False. )
        call Calculate_Vector_Potential( time, aext(:) )
        ! Update the (total) vector potential
        atot(:) = aind(:) + aext(:)
      end if

      ! INDUCED CURRENT
      ! Update the paramagnetic component of the induced current density
      jdia(:) = -atot(:)*chgval/c/omega
      ! Update the paramagnetic component of the induced current density
      jpara(:) = jparanext(:)
      jind(:) = jpara(:)+jdia(:)

      ! HAMILTONIAN
      ham_predcorr(:,:,:) = ham_time(:,:,:)
      call UpdateHam( predcorr=.True., calculateOverlap=.False. )

      ! Check the difference between the two hamiltonians
      err = maxval(abs(ham_predcorr(:,:,:)-ham_time(:,:,:)))
      if ( err <= tolPredCorr ) exit

    end do
    maxStepsReached = (i>maxSteps)
  end subroutine 

  !> Subroutine to initialize all MD related variables
  subroutine init_MD( timeStepRTTDDFT, timeStepMultiplier, molecular_dynamics, &
      MD_outputs, atom_positions, atom_velocities, forces )
    !> Time step used in the real-time TDDFT calculation
    real(dp), intent(in)               :: timeStepRTTDDFT
    !> Integer ratio between the time step used in MD and `timeStepRTTDDFT`
    integer(i32), intent(out)          :: timeStepMultiplier
    !> variable with interfaces to elements defined in the input file
    type(MD_input_keys), intent(inout) :: molecular_dynamics
    !> variable with interfaces to MD outputs
    type(MD_out), intent(inout)        :: MD_outputs
    !> positions of all atoms in cartesian coordinates
    real(dp), allocatable, intent(out) :: atom_positions(:, :) 
    !> velocities of all atoms in cartesian coordinates
    real(dp), allocatable, intent(out) :: atom_velocities(:, :)
    !> forces acting on all atoms
    type(force), intent(out)           :: forces


    timeStepMultiplier = int( molecular_dynamics%time_step/timeStepRTTDDFT )
    molecular_dynamics%time_step = timeStepMultiplier*timeStepRTTDDFT
    
    call MD_allocate_global_arrays( nspecies )
    call MD_evaluate_charge_val
    
    call forces%allocate_arrays( natmtot )
    call force_rttdft( forces, molecular_dynamics%core_corrections, molecular_dynamics%valence_corrections )
    
    allocate( atom_velocities(3, natmtot) )
    call init_atoms_velocities( atom_velocities )

    allocate( atom_positions(3, natmtot) )
    call init_atoms_positions( atom_positions )
    
    if( molecular_dynamics%basis_derivative ) call Update_basis_derivative( atom_velocities, mathcalB, B_time, B_past )
    
    if ( rank == 0 ) then
      call MD_outputs%open_files( natmtot, molecular_dynamics%print_all_force_components  )
      call write_MD_outputs( time, atom_positions, atom_velocities, forces, &
                molecular_dynamics%print_all_force_components, MD_outputs )
    end if

  end subroutine

  !> Initialize the velocities of each atom (needed for Ehrenfest MD)
  subroutine init_atoms_velocities(at_velocities)
    !> Velocities of the nuclei at time \( t = 0 \)
    real(dp), intent(inout)         :: at_velocities(:,:)
    integer(i32) :: is, ia, ias

    call assert( size(at_velocities, 1) == 3, 'at_velocities must have size = 3 along dim = 1' )
    call assert( size(at_velocities, 2) == natmtot, 'at_velocities must have size = natmtot along dim = 2' )

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        at_velocities(1:3,ias) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%velocity(1:3)
      end do
    end do
  end subroutine

  !> Initialize the positions of each atom from the global variable `atposc`. 
  !> Needed for Ehrenfest MD
  subroutine init_atoms_positions(at_positions)
    !> Positions of the nuclei at time \( t = 0 \)
    real(dp), intent(inout)         :: at_positions(:,:)
    integer(i32) :: is, ia, ias

    call assert( size(at_positions, 1) == 3, 'at_positions must have size = 3 along dim = 1' )
    call assert( size(at_positions, 2) == natmtot, 'at_positions must have size = natmtot along dim = 2' )
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        at_positions(1:3, ias) = atposc(1:3, ia, is)
      end do
    end do
  end subroutine

  subroutine write_MD_outputs( t, at_positions, at_velocities, forces, print_all_force_components, MD_outputs )
    !> time \(t\)
    real(dp), intent(in) :: t
    !> Positions of the nuclei at time \( t  \)
    real(dp), intent(in)         :: at_positions(:,:)
    !> Velocities of the nuclei at time \( t \)
    real(dp), intent(in)         :: at_velocities(:,:)
    !> forces acting on all atoms
    type(force), intent(in)      :: forces
    !> if `.True.`, print out all contributions to the total force
    logical, intent(in)          :: print_all_force_components
    !> variable with interfaces to MD outputs
    type(MD_out), intent(inout)  :: MD_outputs
    
    if ( print_all_force_components ) then
      call MD_outputs%write_to_files( t, at_positions, at_velocities, forces%total, &
        forces )
    else 
      call MD_outputs%write_to_files( t, at_positions, at_velocities, forces%total )
    end if
  end subroutine

  subroutine deallocate_global_arrays( deallocate_ehrenfest_arrays )
    !> when `.True.`, also deallocate arrays used in Ehrenfest MD
    logical, intent(in) :: deallocate_ehrenfest_arrays

    if ( allocated(apwalm) ) deallocate( apwalm )
    if ( allocated(evecfv_gnd) ) deallocate( evecfv_gnd )
    if ( allocated(evecsv) ) deallocate( evecsv )
    if ( allocated(evecfv_time) ) deallocate( evecfv_time )
    if ( allocated(overlap) ) deallocate( overlap )
    if ( allocated(ham_time) ) deallocate( ham_time )
    if ( allocated(ham_past) ) deallocate( ham_past )
    if ( allocated(pmat) ) deallocate( pmat )
    if ( predictorCorrector ) deallocate( ham_predcorr, evecfv_save )
    if ( associated(input%xs%realTimeTDDFT%laser) ) then
      if ( nkicks >= 1 ) then
        deallocate( wkick, dirkick, amplkick, t0kick )
      end if
      if ( ntrapcos >= 1 ) then
        deallocate( dirtrapcos, ampltrapcos, omegatrapcos, phasetrapcos )
        deallocate( t0trapcos, trtrapcos, wtrapcos )
      end if
      if ( nsinsq >= 1 ) then
        deallocate( dirsinsq, amplsinsq, omegasinsq )
        deallocate( phasesinsq, t0sinsq, tpulsesinsq )
      end if
    end if
    
    if ( deallocate_ehrenfest_arrays ) then
      call MD_deallocate_global_arrays
      if ( allocated( pmatmt ) ) deallocate( pmatmt )
      if ( allocated( mathcalH ) ) deallocate( mathcalH )
      if ( allocated( mathcalB ) ) deallocate( mathcalB )
      if ( allocated( B_past ) ) deallocate( B_past )
      if ( allocated( B_time ) ) deallocate( B_time )
    end if
  end subroutine
end module rttddft_main
