

! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine initocc(nbf, naf)
  use modmain
use modinput
  use modxs
  implicit none
  ! arguments
  integer :: nbf, naf
  ! local variables
  character(*), parameter :: thisnam='initocc'
  integer :: nvalel
  ! number of valence electrons
  nvalel=nint(chgval/2.d0)
  ! number of states below Fermi energy
  if (nbf.eq.0) nbf=nvalel
  if (nbf.gt.nvalel) then
     write(unitout, '("Warning(", a, "): number of states below Fermi energy too &
     &large - adjusting to number of valence states")') trim(thisnam)
  end if
  ! number of states above Fermi energy
  if (naf.eq.0) then
     ! if "naf" is not specified define it using "nempty"
     naf=input%groundstate%nempty+1
  else
     ! check if number is too large
     if (naf.gt.(input%groundstate%nempty+1)) then
	naf=input%groundstate%nempty+1
	write(unitout, '("Warning(", a, "): number of states above Fermi energy &
	&too large - adjusting using number of empty states")') trim(thisnam)
     end if
  end if
end subroutine initocc
