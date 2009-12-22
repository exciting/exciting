!
!
!
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Logical Function tqgamma (iq)
      Use modmain
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Real (8) :: epsg = 1.d-12
      tqgamma = .False.
      If (sum(Abs(vqc(:, iq))) .Lt. epsg) tqgamma = .True.
End Function tqgamma
