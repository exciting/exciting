
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmwriteengy(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
write(fnum,*)
write(fnum,'("Energies :")')
write(fnum,'(" electronic kinetic",T30,": ",G18.10)') engykn
write(fnum,'(" core electron kinetic",T30,": ",G18.10)') engykncr
write(fnum,'(" Coulomb",T30,": ",G18.10)') engyvcl
write(fnum,'(" Madelung",T30,": ",G18.10)') engymad
write(fnum,'(" exchange-correlation",T30,": ",G18.10)') engyx
if (rdmtemp.gt.0.d0) then
  write(*,'(" entropy",T30,": ",G18.10)') rdmentrpy
end if
write(fnum,'(" total energy",T30,": ",G18.10)') engytot
call flushifc(fnum)
return
end subroutine

