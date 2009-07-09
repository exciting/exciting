


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

integer function idxkkp(ik, ikp, n)
  implicit none
  ! arguments
  integer, intent(in) :: ik, ikp, n
  ! local variables
  if ((ik.le.0).or.(ikp.le.0).or.(n.le.0)) then
     write(*, *)
     write(*, '("Error(idxkkp): negative indices or number of points")')
     write(*, *)
     call terminate
  end if
  if (ik.gt.ikp) then
     write(*, *)
     write(*, '("Error(idxkkp): ik > ikp")')
     write(*, *)
     call terminate
  end if
  ! (i,j) -> (i-1)i/2 + (N-i+1)(i-1) + j-i+1
  idxkkp=-(ik*(ik-1))/2+n*(ik-1)+ikp
end function idxkkp
