subroutine check_acoustic(vec, acoustic)
!
use mod_atoms
use raman_params
!
implicit none
! arguments
real(8) :: vec(3*natmtot)
logical :: acoustic
! local variables
integer :: i
integer :: icomp(3*natmtot)
logical :: acousticx, acousticy, acousticz
!
acoustic = .false.
do i = 1, 3*natmtot
   if (abs(vec(i)) .le. eps) then
      icomp(i) = 0
   elseif (vec(i) .gt. eps) then
      icomp(i) = 1
   else
      icomp(i) = -1
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
