!
!-----------------------------------------------------------------------------80.
! Copyright (C) 2013- exciting team
! This file is distributed under the terms of the GNU General Public License.
! Created       on 15-10-2013 Pasquale Pavone (exciting team)
! Modified      on 15-11-2013 Pasquale Pavone (exciting team)
! Last modified on 14-06-2014 Pasquale Pavone (exciting team)
!-----------------------------------------------------------------------------80
!
subroutine writepositions(fnum,verbosity)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
    integer, intent(in) :: verbosity
! local variables
    integer :: is, ia, i
    character(10) :: acoord

    if (verbosity < 1) return

    i = 0
    acoord = "lattice"
    if (input%structure%cartesian) acoord = "cartesian"
    write(fnum,*)
    write(fnum,'(" Atomic positions at this step (",A,") :")') trim(acoord)
    do is = 1, nspecies
        do ia = 1, natoms (is)
            i = i+1
            if (input%structure%cartesian) then  
                write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
               &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  atposc(:,ia,is)
            else
                write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
               &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            end if
        end do
    end do
    
    return
end subroutine
