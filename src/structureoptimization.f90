!BOP
! !ROUTINE: structureoptimization
! !INTERFACE:
!
!
subroutine structureoptimization
! !USES:
    Use modinput
    Use modmain
    Use modmpi
    Use scl_xml_out_Module
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!   Created February 2013 (DIN)
!EOP
!BOC
    Implicit None
    integer :: is, ia

    if (input%structureoptimization%method=="simple") then
        
        do while (forcemax > input%structureoptimization%epsforce)
            Call updatpos
            If (rank .Eq. 0) Then
                Write (60,*)
                Write (60, '("+--------------------------+")')
                Write (60, '("| Updated atomic positions |")')
                Write (60, '("+--------------------------+")')
                Do is = 1, nspecies
                    Write (60,*)
                    Write (60, '("Species : ", I4, " (", A, ")")') &
                   &  is, trim (input%structure%speciesarray(is)%species%chemicalSymbol)
                    Write (60, '(" atomic positions (lattice) :")')
                    Do ia = 1, natoms (is)
                        Write (60, '(I4, " : ", 3F14.8)') ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                    End Do
                End Do
            End If
! begin new self-consistent loop with updated positions
            call scf_cycle
        end do

    else if (input%structureoptimization%method=="lbfgs") then
        
        If (rank .Eq. 0) Then
            Write (60,*)
            Write (60, '("+---------------------------------------------------------+")')
            Write (60, '("| Use L-BFGS-B method for optimizing the atomic positions |")')
            Write (60, '("+---------------------------------------------------------+")')
        End If
        call lbfgs_driver

    else
      
        write(*,*) 'ERROR(structureoptimization.f90): Unknown method for structure optimization!'
        stop
      
    end if

    return
end subroutine
