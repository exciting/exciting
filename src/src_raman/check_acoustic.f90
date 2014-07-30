! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
subroutine check_acoustic(vec, acoustic)
!
use mod_atoms
use raman_params
!
implicit none
! arguments
complex(8), intent(inout) :: vec(3*natmtot)
logical, intent(out) :: acoustic
! local variables
real(8) :: phi
integer :: i
integer :: icomp(3*natmtot)
logical :: acousticx, acousticy, acousticz
!
acoustic = .false.
do i = 1, 3*natmtot
   ! clean close-to-real numbers
   if (abs(aimag(vec(i))) .lt. eps) vec(i) = cmplx(dble(vec(i)), 0.d0, 8)
   phi = atan2(aimag(vec(i)), dble(vec(i)))
   if (abs(vec(i)) .le. eps) then
      icomp(i) = 0
   else
      if (phi .gt. 0.d0) then
         icomp(i) = 1
      else
         icomp(i) = -1
      endif
   endif
enddo
acousticx = .false.; acousticy = .false.; acousticz = .false.
if (all(icomp(1:3*natmtot-2:3) .eq.  0) .or. &
 &  all(icomp(1:3*natmtot-2:3) .eq.  1) .or. &
 &  all(icomp(1:3*natmtot-2:3) .eq. -1) ) acousticx = .true.
if (all(icomp(2:3*natmtot-1:3) .eq.  0) .or. &
 &  all(icomp(2:3*natmtot-1:3) .eq.  1) .or. &
 &  all(icomp(2:3*natmtot-1:3) .eq. -1) ) acousticy = .true.
if (all(icomp(3:3*natmtot:3) .eq.  0) .or. &
 &  all(icomp(3:3*natmtot:3) .eq.  1) .or. &
 &  all(icomp(3:3*natmtot:3) .eq. -1) ) acousticz = .true.
if (acousticx .and. acousticy .and. acousticz) acoustic = .true.
!
end subroutine
