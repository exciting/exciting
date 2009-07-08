

! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine writefsm(fnum)
use modmain
use modinput
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer::is, ia
if (input%groundstate%spin%fixspinnumber.eq.0) return
write(fnum, *)
if ((input%groundstate%spin%fixspinnumber.eq.1).or.(input%groundstate%spin%fixspinnumber.eq.3)) then
  write(fnum, '("FSM global effective field", T30, ": ", 3G18.10)') bfsmc(1:ndmag)
end if
if ((input%groundstate%spin%fixspinnumber.eq.2).or.(input%groundstate%spin%fixspinnumber.eq.3)) then
  write(fnum, '("FSM local muffin-tin effective fields :")')
  do is=1, nspecies
    write(fnum, '(" species : ", I4, " (", A, ")")') is, trim(input%structure%speciesarray(is)%species%chemicalSymbol)
    do ia=1, natoms(is)
      write(fnum, '("  atom ", I4, T30, ": ", 3G18.10)') ia, bfsmcmt(1:ndmag, ia, is)
    end do
  end do
end if
return
end subroutine
