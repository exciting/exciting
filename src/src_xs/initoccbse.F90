
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initoccbse(nbf,naf)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer :: nbf,naf
  ! local variables
  character(*), parameter :: thisnam='initoccbse'
  integer :: nvalel
  ! number of valence electrons
  nvalel=nint(chgval/2.d0)
  ! number of states below Fermi energy
  if (nbf.eq.-1) nbf=nvalel
  if (nbf.gt.nvalel) then
     write(unitout,'(a,2i6)') 'Error('//trim(thisnam)//'): number of &
          &states below Fermi energy for BSE too large (proposed/max):',&
          nbf,nvalel
     call terminate
  end if
  ! number of states above Fermi energy
  if (naf.eq.-1) then
     ! if "naf" is not specified define it using "nempty"
     naf=nempty+1
  else
     ! check if number is too large
     if (naf.gt.(nempty+1)) then
        naf=nempty+1
        write(unitout,'("Info(",a,"): naf too large: adjusting &
             & to ",I6)') naf
     end if
  end if
end subroutine initoccbse
