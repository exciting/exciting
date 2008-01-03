
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getridx(procs,set,i,ir)
  ! Get index ir relative to partition determined by overall index i.
  use modmpi, only: firstofset,procofindex
  implicit none
  ! arguments
  integer :: procs
  integer, intent(in) :: set,i
  integer, intent(out) :: ir
  ! local variables
  integer :: procst
  ! executable statements begin
  procst=procs
  ir=i-firstofset(procofindex(i,set),set)+1
end subroutine getridx
