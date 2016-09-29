! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: getridx
! !INTERFACE:
subroutine getridx(set, i, ir)
! !USES:
  use modmpi, only: firstofset, procofindex
! !INPUT/OUTPUT PARAMETERS:
! IN:
! integer :: set ! Number of elements distributed over MPI
! integer :: i   ! Absolut element index
! OUT:
! integer :: ir  ! Relative element index with respect to the set of elements
!                ! the current rank was given
! !DESCRIPTION:
!   Get index ir relative to partition determined by overall index i.
! !REVISION HISTORY:
!   Added to documentation scheme and removed dead code parts. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: set, i
  integer, intent(out) :: ir

  ! Executable statements begin
  ir = i - firstofset(procofindex(i, set), set) + 1
end subroutine getridx
!EOC
