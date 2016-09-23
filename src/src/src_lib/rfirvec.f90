!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
Real (8) Function rfirvec (ngrid, ainv, vc, rfir)
      Implicit None
! arguments
      Integer, Intent (In) :: ngrid (3)
      Real (8), Intent (In) :: ainv (3, 3)
      Real (8), Intent (In) :: vc (3)
      Real (8), Intent (In) :: rfir (ngrid(1), ngrid(2), ngrid(3))
! local variables
      Integer :: i1, i2, i3, j1, j2, j3
      Real (8), Parameter :: eps = 1.d-6
      Real (8) :: v1, v2, v3, p1, p2, p3, q1, q2, q3
      Real (8) :: f00, f01, f10, f11, f0, f1
! input vector in lattice coordinates
      v1 = ainv (1, 1) * vc (1) + ainv (1, 2) * vc (2) + ainv (1, 3) * &
     & vc (3)
      v2 = ainv (2, 1) * vc (1) + ainv (2, 2) * vc (2) + ainv (2, 3) * &
     & vc (3)
      v3 = ainv (3, 1) * vc (1) + ainv (3, 2) * vc (2) + ainv (3, 3) * &
     & vc (3)
! map lattice coordinates to [0,1) interval
      i1 = Int (v1)
      i2 = Int (v2)
      i3 = Int (v3)
      v1 = v1 - dble (i1)
      v2 = v2 - dble (i2)
      v3 = v3 - dble (i3)
      If (v1 .Lt. 0.d0) v1 = v1 + 1.d0
      If (v2 .Lt. 0.d0) v2 = v2 + 1.d0
      If (v3 .Lt. 0.d0) v3 = v3 + 1.d0
      If (1.d0-v1 .Lt. eps) v1 = 0.d0
      If (1.d0-v2 .Lt. eps) v2 = 0.d0
      If (1.d0-v3 .Lt. eps) v3 = 0.d0
      If (v1 .Lt. eps) v1 = 0.d0
      If (v2 .Lt. eps) v2 = 0.d0
      If (v3 .Lt. eps) v3 = 0.d0
! determine coordinates on grid
      v1 = dble (ngrid(1)) * v1
      v2 = dble (ngrid(2)) * v2
      v3 = dble (ngrid(3)) * v3
      i1 = Int (v1)
      i2 = Int (v2)
      i3 = Int (v3)
! use trilinear interpolation with neighbouring points
      p1 = v1 - dble (i1)
      p2 = v2 - dble (i2)
      p3 = v3 - dble (i3)
      q1 = 1.d0 - p1
      q2 = 1.d0 - p2
      q3 = 1.d0 - p3
      i1 = modulo (i1, ngrid(1)) + 1
      i2 = modulo (i2, ngrid(2)) + 1
      i3 = modulo (i3, ngrid(3)) + 1
      j1 = modulo (i1, ngrid(1)) + 1
      j2 = modulo (i2, ngrid(2)) + 1
      j3 = modulo (i3, ngrid(3)) + 1
      f00 = rfir (i1, i2, i3) * q1 + rfir (j1, i2, i3) * p1
      f01 = rfir (i1, i2, j3) * q1 + rfir (j1, i2, j3) * p1
      f10 = rfir (i1, j2, i3) * q1 + rfir (j1, j2, i3) * p1
      f11 = rfir (i1, j2, j3) * q1 + rfir (j1, j2, j3) * p1
      f0 = f00 * q2 + f10 * p2
      f1 = f01 * q2 + f11 * p2
      rfirvec = f0 * q3 + f1 * p3
      Return
End Function
