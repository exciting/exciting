! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: screen
! !INTERFACE:
subroutine screen
! !USES:
  use modmpi
  use modxs, only: nwdf, unitout
  use m_genfilname
  use m_filedel
! !DESCRIPTION:
! This is a wrapper subroutine that backs up the number of 
! $ \omega $ points and calls the {\tt df} subroutine. Afterwards
! it restores them. It also sets the output file extensions for
! screening related quantities.
! 
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! local variables
  integer :: nwdft

  ! Back up number of energy points
  nwdft = nwdf

  ! Change file extension variable in mod_misc to '_SCR.OUT' for the
  ! following calculations
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Call dielectric function with only one frequency point
  ! df is a wrapper for dfq(iq) which for 'screen' sets nwdf=1, i.e. static screening
  call df

  ! Alternative for checking only:
  nwdf = nwdft

  if(rank == 0) then
    write(unitout, '(a)') "Info(screen): Screening finished"
  end if
end subroutine screen
!EOC
