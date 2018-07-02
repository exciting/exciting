
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine groundstatetasklauncher
    Use modinput
    Use modmain, Only: task,xctype
    Use inputdom
    use modmpi
    Implicit None
    
    call delete_warnings
    splittfile= .true.
    If ( .Not. (associated(input%groundstate%solver))) Then
        ! set the default values if tddft element not present
        input%groundstate%solver => getstructsolver (emptynode)
    End If
    If (((xctype(2) .Ge. 400).Or. (xctype(1) .Ge. 400)).and. &
    & (.not.associated(input%groundstate%Hybrid))) Then 
           input%groundstate%Hybrid => getstructHybrid (emptynode)
    End If  
    If  (associated(input%groundstate%Hybrid)) Then
        If (input%groundstate%Hybrid%exchangetypenumber .Eq. 2) Then
            If (.not.associated(input%groundstate%OEP)) Then
               input%groundstate%OEP => getstructOEP (emptynode)
            End If 
        End If                
    End If
    If (input%groundstate%xctypenumber .Lt. 0) Then
        If (.not.associated(input%groundstate%OEP)) Then
           input%groundstate%OEP => getstructOEP (emptynode)
        End If 
    End If            
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
        ! DFT / OEP
        Else If (associated(input%groundstate%Hybrid)) Then
                If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
                    Call hybrids
                Else 
                    Call gndstate
                End If    
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

end subroutine
