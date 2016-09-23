!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: vecfbz
! !INTERFACE:
!
!
Subroutine vecfbz (eps, bvec, vpl, iv)
! !INPUT/OUTPUT PARAMETERS:
!   eps  : zero component tolerance (in,real)
!   bvec : reciprocal lattice vectors (in,real(3,3))
!   vpl  : input vector in lattice coordinates (inout,real(3))
!   iv   : integer parts of vpl (out,integer(3))
! !DESCRIPTION:
!   Maps a vector in lattice coordinates to the first Brillouin zone. This is
!   done by first removing its integer components and then adding primitive
!   reciprocal lattice vectors until the shortest vector is found.
!
! !REVISION HISTORY:
!   Created September 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: eps
      Real (8), Intent (In) :: bvec (3, 3)
      Real (8), Intent (Inout) :: vpl (3)
      Integer, Intent (Out) :: iv (3)
! local variables
      Integer :: i1, i2, i3, j1, j2, j3
      Real (8) :: v0 (3), v1 (3), v2 (3), v3 (3), t1, t2
! map vector to [0,1) interval
      Call r3frac (eps, vpl, iv)
      v0 (:) = bvec (:, 1) * vpl (1) + bvec (:, 2) * vpl (2) + bvec (:, &
     & 3) * vpl (3)
      t1 = v0 (1) ** 2 + v0 (2) ** 2 + v0 (3) ** 2
      j1 = 0
      j2 = 0
      j3 = 0
      Do i1 = - 1, 0
         v1 (:) = v0 (:) + dble (i1) * bvec (:, 1)
         Do i2 = - 1, 0
            v2 (:) = v1 (:) + dble (i2) * bvec (:, 2)
            Do i3 = - 1, 0
               v3 (:) = v2 (:) + dble (i3) * bvec (:, 3)
               t2 = v3 (1) ** 2 + v3 (2) ** 2 + v3 (3) ** 2
               If (t2 .Lt. t1+eps) Then
                  j1 = i1
                  j2 = i2
                  j3 = i3
                  t1 = t2
               End If
            End Do
         End Do
      End Do
      vpl (1) = vpl (1) + dble (j1)
      vpl (2) = vpl (2) + dble (j2)
      vpl (3) = vpl (3) + dble (j3)
      iv (1) = iv (1) - j1
      iv (2) = iv (2) - j2
      iv (3) = iv (3) - j3
      Return
End Subroutine
!EOC
