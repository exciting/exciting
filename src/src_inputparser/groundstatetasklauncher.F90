
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
module groundstatetasklauncher
  use modinput, only: input
  use precision, only: dp, i32

  implicit none
  private

  public :: launch_groundstate

contains
subroutine launch_groundstate
    use cdft, only: cdft_gs_run_request, file_extension_CDFT, file_extension_GS, gs_run_before_CDFT, &
        set_status_to_finished_CDFT, set_status_to_running_CDFT, single_shot, skip
    use dfthalf, only: run_dft_half_nscf
    use inputdom, only: emptynode
    use mod_misc, only: filext, task
    use mod_selfconsistent_gw, only: prepare_current_iteration, check_convergence, iteration, selfconsistent_gw_eps, &
                                     initialize_selfconsitent_gw, gw_first_iteration
    use mod_potential_and_density, only: xctype
    Use modinput, only: getstructHybrid, getstructOEP, getstructsolver
    use modmpi, only: terminate_if_false, splittfile
    

    integer(i32) :: task_backup, maxscl_backup
    integer(kind( gs_run_before_CDFT )) :: cdft_gs_run
    character(len=:), allocatable :: string
    logical :: is_cdft_calculation, is_hybrid_calculation, is_exchange_hartree_fock
    logical :: dft_half_nscf

    !
    ! If this GS run is within an external QSGW cycle 
    ! this calls manage them
    !
    ! 1. Prepare all the information
    call initialize_selfconsitent_gw()
    ! 2. Execute the GS with the proper 
    !    globals. This only happens beyond the
    !    first iteration. For the latter
    !    the calculation is deferred to 
    !    the usual.
    call prepare_current_iteration()
    ! 3. Check the convergence of the outer cycle
    call check_convergence(iteration, selfconsistent_gw_eps)

    call delete_warnings
    splittfile= .true.
    If ( .Not. (associated(input%groundstate%solver))) Then
        ! set the default values if tddft element not present
        input%groundstate%solver => getstructsolver (emptynode)
    End If
    ! default (HF-based) hybrid functionals
    If ((xctype(1) >= 400).and.(.not.associated(input%groundstate%Hybrid))) Then
        input%groundstate%Hybrid => getstructHybrid(emptynode)
    End If
    ! EXX-OEP-based hybrid functionals
    If (associated(input%groundstate%Hybrid)) Then
        If (input%groundstate%Hybrid%exchangetypenumber .Eq. 2) Then
            If (.not.associated(input%groundstate%OEP)) Then
               input%groundstate%OEP => getstructOEP (emptynode)
            End If
        End If
    End If
    If (input%groundstate%xctypenumber .Lt. 0) Then
        If (.not.associated(input%groundstate%OEP)) Then
           input%groundstate%OEP => getstructOEP(emptynode)
        End If
    End If

    ! Interface to input elements defined for a hybrid calculation
    is_hybrid_calculation = associated(input%groundstate%Hybrid)
    is_exchange_hartree_fock = .false.
    if ( is_hybrid_calculation ) is_exchange_hartree_fock = ( input%groundstate%Hybrid%exchangetypenumber == 1 )

    ! Interface to input elements defined for a constrained DFT calculation
    is_cdft_calculation = associated(input%groundstate%constrainedDFT)
    if ( is_cdft_calculation ) cdft_gs_run = cdft_gs_run_request( input%groundstate%constrainedDFT%groundstateRun )

    ! Interface to input elements for DFT-1/2
    dft_half_nscf = associated( input%groundstate%dfthalf )
    if( dft_half_nscf ) dft_half_nscf = input%groundstate%dfthalf%NSCF

    If (input%groundstate%do .Eq. "fromscratch") Then
        If (associated(input%relax)) Then
            task = 2
        Else
            task = 0
        End If
    Else
        If (associated(input%relax)) Then
            task = 3
        Else
            task = 1
        End If
    End If
    If (associated(input%groundstate%OEP)) Then
        input%groundstate%scfconv = 'potential'
    End If
    If (input%groundstate%do .Ne. "skip") then
        ! Hartree Fock
        If  (associated(input%groundstate%HartreeFock)) Then
            task = 5
            Call hartfock

        ! Constrained DFT calculation
        Else If ( is_cdft_calculation ) then
            string = filext
            if ( is_hybrid_calculation ) &
                call terminate_if_false( is_exchange_hartree_fock , "The combination of OEP and CDFT is not implemented yet." )
            if ( cdft_gs_run /= skip ) then
                filext = file_extension_GS
                if ( cdft_gs_run == single_shot ) then
                    task_backup = task
                    maxscl_backup = input%groundstate%maxscl
                    task = 1 ! start from existing `STATE.OUT`
                    input%groundstate%maxscl = 1 ! single-shot run
                end if
                call hybrid_or_gndstate( is_hybrid_calculation .and. is_exchange_hartree_fock )
                if ( cdft_gs_run == single_shot ) then
                    task = task_backup
                    input%groundstate%maxscl = maxscl_backup
                end if
            end if
            if ( task == 7 ) task = 0 ! for hybrid calculation
            filext = file_extension_CDFT
            call set_status_to_running_CDFT()
            call hybrid_or_gndstate( is_hybrid_calculation .and. is_exchange_hartree_fock )
            call set_status_to_finished_CDFT()
            filext = string

        ! DFT / OEP
        Else If (is_hybrid_calculation) Then
            call hybrid_or_gndstate(is_exchange_hartree_fock)
        Else
            Call gndstate
        End If
        ! do conversion to XML format if requested
        if (associated(input%groundstate%output)) then
             if (input%groundstate%output%state .eq. "XML") call portstate(1)
        end if

    else
        splittfile= .False.
    end if
    
    if( dft_half_nscf ) then
        call run_dft_half_nscf
    end if
end subroutine

subroutine hybrid_or_gndstate(is_hartree_fock)
    logical, intent(in) :: is_hartree_fock
    if (is_hartree_fock)  then
        call hybrids
    else
        call gndstate
    end if
end subroutine hybrid_or_gndstate



end module