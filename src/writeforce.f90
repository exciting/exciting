!
! Copyright (C) 2004-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine writeforce(fnum)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
! local variables
    integer :: is, ia, ias
    real(8) :: t1
    write(fnum,'("Atomic forces:")')
    do is = 1, nspecies
        do ia = 1, natoms (is)
            ias = idxas (ia, is)
            write(fnum,'(" atom ",I4,4x,A2,T30)') &
           &  ia, trim(input%structure%speciesarray(is)%species%chemicalSymbol)
            write(fnum,'(4x,"total force",T30,": ",3F14.8)') forcetot(:,ias)
            if ( (input%groundstate%outputlevelnumber>1).or. &
                 (input%structureoptimization%outputlevelnumber>1) ) then
                write(fnum,'(4x,"Hellmann-Feynman",T30,": ",3F14.8)') forcehf(:,ias)
                write(fnum,'(4x,"core correction",T30,": ",3F14.8)') forcecr(:,ias)
                write(fnum,'(4x,"IBS",T30,": ",3F14.8)') forceibs(:,ias)
                t1 = Sqrt (forcetot(1, ias)**2+forcetot(2,ias)**2+forcetot(3, ias)**2)
                write(fnum,'(4x,"total magnitude",T30,": ",F14.8)') t1
            end if
        end do
    end do
    return
end subroutine
