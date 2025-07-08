
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine groundstatetasklauncher
    use cdft, only: file_extension_CDFT, file_extension_GS, set_status_to_finished_CDFT, set_status_to_running_CDFT
    Use modinput
    Use modmain,     only: task, xctype
    Use inputdom
    use modmpi, only: terminate_if_false, splittfile
    use mod_misc, only: filext

    Implicit None

    character(len=:), allocatable :: string
    logical :: is_cdft_calculation, skip_gnd_in_cdft_calculation
    logical :: is_hybrid_calculation, is_exchange_hartree_fock

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
    skip_gnd_in_cdft_calculation = .false.
    if ( is_cdft_calculation ) skip_gnd_in_cdft_calculation = input%groundstate%constrainedDFT%skipgnd

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
        Else If ( is_cdft_calculation ) Then
            string = filext
            if ( is_hybrid_calculation ) then
                call terminate_if_false(  is_exchange_hartree_fock , " ERROR: The combination of OEP and CDFT is not implemented yet." )
            end if
            if ( .not. skip_gnd_in_cdft_calculation ) then
                filext = file_extension_GS
                call hybrid_or_gndstate( is_hybrid_calculation .and. is_exchange_hartree_fock )
            end if
            if (task==7) then  ! for hybrid calculation
                task=0
            end if
            filext = file_extension_CDFT
            call set_status_to_running_CDFT
            call hybrid_or_gndstate( is_hybrid_calculation .and. is_exchange_hartree_fock )
            call set_status_to_finished_CDFT
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

contains
    subroutine hybrid_or_gndstate(is_hartree_fock)
        logical, intent(in) :: is_hartree_fock
        if (is_hartree_fock)  then
            call hybrids
        else
            call gndstate
        end if
    end subroutine hybrid_or_gndstate

end subroutine
