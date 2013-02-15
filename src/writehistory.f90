!
!BOP
! !ROUTINE: writehistory
! !INTERFACE:
!
Subroutine writehistory
! !USES:
    Use modinput
    Use modmain
!
! !DESCRIPTION:
!   Outputs the atomic positions and forces to file, in the XYZ format,
!   which may be then used to visualize the optimization trajectory
!
! !REVISION HISTORY:
!   
!EOP
!BOC
    Implicit None
    Integer :: is, ia, ias
    Real (8) :: v (3)
    character(80) :: historyformat
    
    if (input%structureoptimization%history) then

        historyformat='xyz'
        select case (historyformat)
        case('xyz','XYZ')
            Open (51, File='history.xyz', Action='WRITE', POSITION='APPEND')
            Write (51,*) natmtot
            Write (51,*) "Total Energy =",  engytot*27.211396641344194, " eV "
            Do is = 1, nspecies
                Do ia = 1, natoms (is)
                    ias = idxas (ia, is)
                    v(:) = atposc(:,ia,is)
                    Write (51,'(A6, 3F14.8, 3F14.8)')(input%structure%speciesarray(is)%species%chemicalSymbol), &
                    &    v(:)*0.529177249, forcetot(:, ias)*51.42208245
                End Do
            End Do
            Close (51)
        case default
            write(*,*)'ERROR(writehistory.f90): Unknown output format'
            stop
        end select
    
    end if

    Return
End Subroutine
!EOC
