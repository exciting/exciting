!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine potnucl (ptnucl, nr, r, zn, vn)
      Implicit None
! arguments
      Logical, Intent (In) :: ptnucl
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: zn
      Real (8), Intent (Out) :: vn (nr)
! local variables
      Integer :: ir
! nuclear radius constant in Bohr
      Real (8), Parameter :: r0 = 1.25d-15 / 0.52917720859d-10
      Real (8) :: rn, t1, t2
      If (zn .Eq. 0.d0) Then
         vn (:) = 0.d0
         Return
      End If
      If (ptnucl) Then
! nucleus is taken to be a point particle
         Do ir = 1, nr
            vn (ir) = zn / r (ir)
         End Do
      Else
! nucleus has a finite radius approximated by r0*A^(1/3)
         rn = r0 * Abs (zn) ** (1.d0/3.d0)
         t1 = zn / (2.d0*rn**3)
         t2 = 3.d0 * rn ** 2
         Do ir = 1, nr
            If (r(ir) .Lt. rn) Then
               vn (ir) = t1 * (t2-r(ir)**2)
            Else
               vn (ir) = zn / r (ir)
            End If
         End Do
      End If
      Return
End Subroutine
