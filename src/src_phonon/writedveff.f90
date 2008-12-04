
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writedveff(iq,is,ia,ip,dveffmt,dveffir)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ip
complex(8), intent(in) :: dveffmt(lmmaxvr,nrcmtmax,natmtot0)
complex(8), intent(in) :: dveffir(ngrtot0)
! local variables
integer js
character(256) fext
call phfext(iq,is,ia,ip,fext)
open(50,file='DVEFF'//trim(fext),action='WRITE',form='UNFORMATTED')
write(50) version
write(50) nspecies
write(50) lmmaxvr
do js=1,nspecies
  write(50) natoms0(js)
  write(50) nrcmt(js)
end do
write(50) ngrid0
write(50) dveffmt,dveffir
close(50)
return
end subroutine

