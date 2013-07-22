!
! Copyright (C) 2004-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine writeforce(fnum,verbosity)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
    integer, intent(in) :: verbosity
! local variables
    integer :: is, ia, ias
    real(8) :: t1
    
    if (verbosity > 1) then
        if (input%groundstate%tfibs) then
            write(fnum,'("Atomic force components (including IBS):")')
        else
            write(fnum,'("Atomic force components:")')
        end if
    else
        if (input%groundstate%tfibs) then
            write(fnum,'("Total atomic forces:")')
        else
            write(fnum,'("Total atomic forces (including IBS):")')
        end if
    end if
    do is = 1, nspecies
        do ia = 1, natoms (is)
            ias = idxas (ia, is)
            if (verbosity > 1) then
                write(fnum,'("  atom ",I4,4x,A2)') &
               &  ia, trim(input%structure%speciesarray(is)%species%chemicalSymbol)
                write(fnum,'("    Hellmann-Feynman",T30,": ",3F14.8)') forcehf(:,ias)
                write(fnum,'("    core correction",T30,": ",3F14.8)') forcecr(:,ias)
                if (input%groundstate%tfibs) write(fnum,'("    IBS",T30,": ",3F14.8)') forceibs(:,ias)
                t1 = Sqrt (forcetot(1, ias)**2+forcetot(2,ias)**2+forcetot(3, ias)**2)
                write(fnum,'("    total magnitude",T30,": ",F14.8)') t1
            else
                write(fnum,'("  atom ",I4,4x,A2,T30": ",3F14.8)') &
               &  ia, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  forcetot(:,ias)
            end if
        end do
    end do
    return
end subroutine
