!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rotaxang
! !INTERFACE:
!
!
Subroutine rotaxang (eps, rot, det, v, th)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero vector tolerance (in,real)
!   rot : rotation matrix (in,real(3,3))
!   det : matrix determinant (out,real)
!   v   : normalised axis vector (out,real(3))
!   th  : rotation angle (out,real)
! !DESCRIPTION:
!   Given a rotation matrix
!   $$ R(\hat{\bf v},\theta)=
!     \left(\begin{matrix}
!     \cos\theta+x^2(1-\cos\theta) &
!     xy(1-\cos\theta)+z\sin\theta &
!     xz(1-\cos\theta)-y\sin\theta \\
!     xy(1-\cos\theta)-z\sin\theta &
!     \cos\theta+y^2(1-\cos\theta) &
!     yz(1-\cos\theta)+x\sin\theta \\
!     xz(1-\cos\theta)+y\sin\theta &
!     yz(1-\cos\theta)-x\sin\theta &
!     \cos\theta+z^2(1-\cos\theta)
!    \end{matrix}\right), $$
!   this routine determines the axis of rotation $\hat{\bf v}$ and the angle of
!   rotation $\theta$. If $R$ corresponds to an improper rotation then only the
!   proper part is used and {\tt det} is set to $-1$.
!
! !REVISION HISTORY:
!   Created Decmeber 2006 (JKD)
!   Changed "intent(inout)" to "intent(in)" for argument "rot", 2009 (Sagmeister)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: eps
      Real (8), Intent (In) :: rot (3, 3)
      Real (8), Intent (Out) :: det
      Real (8), Intent (Out) :: v (3)
      Real (8), Intent (Out) :: th
! local variables
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8) :: p (3, 3), t1, t2
! external functions
      Real (8) :: r3mdet
      External r3mdet
! find the determinant
      det = rot (1, 2) * rot (2, 3) * rot (3, 1) - rot (1, 3) * rot (2, &
     & 2) * rot (3, 1) + rot (1, 3) * rot (2, 1) * rot (3, 2) - rot (1, &
     & 1) * rot (2, 3) * rot (3, 2) + rot (1, 1) * rot (2, 2) * rot (3, &
     & 3) - rot (1, 2) * rot (2, 1) * rot (3, 3)
      If (Abs(det-1.d0) .Lt. eps) Then
         det = 1.d0
      Else If (Abs(det+1.d0) .Lt. eps) Then
         det = - 1.d0
      Else
         Go To 10
      End If
! proper rotation matrix
      p (:, :) = det * rot (:, :)
      v (1) = (p(2, 3)-p(3, 2)) / 2.d0
      v (2) = (p(3, 1)-p(1, 3)) / 2.d0
      v (3) = (p(1, 2)-p(2, 1)) / 2.d0
      t1 = Sqrt (v(1)**2+v(2)**2+v(3)**2)
      t2 = (p(1, 1)+p(2, 2)+p(3, 3)-1.d0) / 2.d0
      If (Abs(Abs(t2)-1.d0) .Gt. eps) Then
! theta not equal to 0 or pi
         th = Atan2 (t1, t2)
         v (:) = v (:) / t1
      Else
! special case of sin(th)=0
         If (t2 .Gt. 0.d0) Then
! zero angle: axis arbitrary
            th = 0.d0
            v (:) = 1.d0 / Sqrt (3.d0)
         Else
! rotation by pi
            th = pi
            If ((p(1, 1) .Ge. p(2, 2)) .And. (p(1, 1) .Ge. p(3, 3))) &
           & Then
               If (p(1, 1) .Lt. (-1.d0+eps)) Go To 10
               v (1) = Sqrt (Abs(p(1, 1)+1.d0)/2.d0)
               v (2) = (p(2, 1)+p(1, 2)) / (4.d0*v(1))
               v (3) = (p(3, 1)+p(1, 3)) / (4.d0*v(1))
            Else If ((p(2, 2) .Ge. p(1, 1)) .And. (p(2, 2) .Ge. p(3, &
           & 3))) Then
               If (p(2, 2) .Lt. (-1.d0+eps)) Go To 10
               v (2) = Sqrt (Abs(p(2, 2)+1.d0)/2.d0)
               v (3) = (p(3, 2)+p(2, 3)) / (4.d0*v(2))
               v (1) = (p(1, 2)+p(2, 1)) / (4.d0*v(2))
            Else
               If (p(3, 3) .Lt. (-1.d0+eps)) Go To 10
               v (3) = Sqrt (Abs(p(3, 3)+1.d0)/2.d0)
               v (1) = (p(1, 3)+p(3, 1)) / (4.d0*v(3))
               v (2) = (p(2, 3)+p(3, 2)) / (4.d0*v(3))
            End If
         End If
      End If
      Return
10    Continue
      Write (*,*)
      Write (*, '("Error(rotaxang): invalid rotation matrix:")')
      Write (*, '(3G18.10)') rot (1, :)
      Write (*, '(3G18.10)') rot (2, :)
      Write (*, '(3G18.10)') rot (3, :)
      Write (*,*)
      Stop
End Subroutine
!EOC
