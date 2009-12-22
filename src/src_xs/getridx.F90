!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getridx (procs, set, i, ir)
  ! Get index ir relative to partition determined by overall index i.
      Use modmpi, Only: firstofset, procofindex
      Implicit None
  ! arguments
      Integer :: procs
      Integer, Intent (In) :: set, i
      Integer, Intent (Out) :: ir
  ! local variables
      Integer :: procst
  ! executable statements begin
      procst = procs
      ir = i - firstofset (procofindex(i, set), set) + 1
End Subroutine getridx
