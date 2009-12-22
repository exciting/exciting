!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: axangsu2
!
!
Subroutine axangsu2 (v, th, su2)
! !INPUT/OUTPUT PARAMETERS:
!   v   : rotation axis vector (in,real(3))
!   th  : rotation angle (in,real)
!   su2 : SU(2) representation of rotation (out,complex(2,2))
! !DESCRIPTION:
!   Finds the complex ${\rm SU}(2)$ representation of a rotation defined by an
!   axis vector $\hat{\bf v}$ and angle $\theta$. The spinor rotation matrix is
!   given explicitly by
!   $$ R^{1/2}(\hat{\bf v},\theta)=I\cos\frac{\theta}{2}
!    +i(\hat{\bf v}\cdot\vec{\sigma})\sin\frac{\theta}{2}. $$
!
! !REVISION HISTORY:
!   Created August 2007 (JKD)
!   Reversed rotation direction, February 2008 (L. Nordstrom)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: v (3)
      Real (8), Intent (In) :: th
      Complex (8), Intent (Out) :: su2 (2, 2)
! local variables
      Real (8), Parameter :: eps = 1.d-6
      Real (8) :: vn (3), cs, sn, t1
      t1 = Sqrt (v(1)**2+v(2)**2+v(3)**2)
      If (t1 .Lt. eps) Then
         Write (*,*)
         Write (*, '("Error(axangsu2): zero length axis vector")')
         Write (*,*)
         Stop
      End If
! normalise the vector
      vn (:) = v (:) / t1
      cs = Cos (th/2.d0)
      sn = Sin (th/2.d0)
      su2 (1, 1) = cmplx (cs, vn(3)*sn, 8)
      su2 (1, 2) = cmplx (vn(2)*sn, vn(1)*sn, 8)
      su2 (2, 1) = cmplx (-vn(2)*sn, vn(1)*sn, 8)
      su2 (2, 2) = cmplx (cs,-vn(3)*sn, 8)
      Return
End Subroutine
!EOC
