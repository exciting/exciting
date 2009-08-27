


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gentetlinkp
! !INTERFACE:


subroutine gentetlinkp(vpl, tqw)
  ! !USES:
use modinput
  use modmain
  use modtetra
! !DESCRIPTION:
!  Interface to the {\tt gentetlink} routine.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3)
  integer, intent(in) :: tqw
  ! call to interface routine
  call gentetlink(vpl, tqw, input%structure%epslat, bvec, input%groundstate%ngridk, input%groundstate%vkloff, nkpt, &
    &nkptnr, vklnr, &
       ikmapnr)
end subroutine gentetlinkp
!EOC
