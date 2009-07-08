

! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

logical function tqgamma(iq)
  use modmain
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  real(8) :: epsg=1.d-12
  tqgamma=.false.
  if (sum(abs(vqc(:, iq))).lt.epsg) tqgamma=.true.
end function tqgamma
