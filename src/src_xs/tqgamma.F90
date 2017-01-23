! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
logical function tqgamma(iq)
  use mod_qpoint, only: vqc

  ! Arguments
  integer, intent(in) :: iq

  ! Local variables
  real(8), parameter :: epsg = 1.d-12

  tqgamma = .false.
  if(sum(abs(vqc(:, iq))) .lt. epsg) tqgamma = .true.

end function tqgamma
