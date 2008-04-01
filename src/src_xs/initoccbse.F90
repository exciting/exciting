
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initoccbse
  use modmain
  use modxs
  implicit none
  ! local variables
  character(*), parameter :: thisnam='initoccbse'
  integer :: nemptyt, nvalel
  ! number of valence electrons
  nvalel=nint(chgval/2.d0)
  ! number of states below Fermi energy
  if (nstbef.eq.-1) nstbef=nvalel
  if (nstbef.gt.nvalel) then
     write(unitout,'(a,2i6)') 'Error('//trim(thisnam)//'): number of &
          &states below Fermi energy for BSE too large (proposed/max):',&
          nstbef,nvalel
     call terminate
  end if
  ! number of states above Fermi energy
  if (nstabf.eq.-1) then
     ! if "nstabf" is not specified define it using "nempty"
     nstabf=nempty+1
  else
     ! check if number is too large
     if (nstabf.gt.(nempty+1)) then
        nstabf=nempty+1
        write(unitout,'("Info(",a,"): nstabf too large: adjusting &
             &" to ",I6)') nstabf
     end if
  end if
end subroutine initoccbse
