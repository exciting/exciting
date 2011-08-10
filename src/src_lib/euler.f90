!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: euler
! !INTERFACE:
!
!
Subroutine euler (rot, ang)
! !INPUT/OUTPUT PARAMETERS:
!   rot : rotation matrix (in,real(3,3))
!   ang : euler angles (alpha, beta, gamma) (out,real(3))
! !DESCRIPTION:
!   Given a rotation matrix
!   \begin{align*}
!    &R(\alpha,\beta,\gamma)=\\
!    &\left(\begin{matrix}
!     \cos\gamma\cos\beta\cos\alpha-\sin\gamma\sin\alpha &
!     \cos\gamma\cos\beta\sin\alpha+\sin\gamma\cos\alpha &
!    -\cos\gamma\sin\beta \\
!    -\sin\gamma\cos\beta\cos\alpha-\cos\gamma\sin\alpha &
!    -\sin\gamma\cos\beta\sin\alpha+\cos\gamma\cos\alpha &
!     \sin\gamma\sin\beta \\
!     \sin\beta\cos\alpha &
!     \sin\beta\sin\alpha &
!     \cos\beta
!    \end{matrix}\right),
!   \end{align*}
!   this routine determines the Euler angles, $(\alpha,\beta,\gamma)$. This
!   corresponds to the so-called ``y-convention'', which involves the following
!   successive rotations of the coordinate system:
!   \begin{itemize}
!    \item[1.]{The $x_1'$-, $x_2'$-, $x_3'$-axes are rotated anticlockwise
!     through an angle $\alpha$ about the $x_3$ axis}
!    \item[2.]{The $x_1''$-, $x_2''$-, $x_3''$-axes are rotated anticlockwise
!     through an angle $\beta$ about the $x_2'$ axis}
!    \item[3.]{The $x_1'''$-, $x_2'''$-, $x_3'''$-axes are rotated anticlockwise
!     through an angle $\gamma$ about the $x_3''$ axis}
!   \end{itemize}
!   Note that the Euler angles are not necessarily unique for a given rotation
!   matrix.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: rot (3, 3)
      Real (8), Intent (Out) :: ang (3)
! local variables
      Real (8), Parameter :: eps = 1.d-7
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: twopi = 6.2831853071795864769d0
      Real (8) :: det
! find the determinant
      det = rot (1, 2) * rot (2, 3) * rot (3, 1) - rot (1, 3) * rot (2, &
     & 2) * rot (3, 1) + rot (1, 3) * rot (2, 1) * rot (3, 2) - rot (1, &
     & 1) * rot (2, 3) * rot (3, 2) + rot (1, 1) * rot (2, 2) * rot (3, &
     & 3) - rot (1, 2) * rot (2, 1) * rot (3, 3)
      If ((det .Lt. 1.d0-eps) .Or. (det .Gt. 1.d0+eps)) Then
         Write (*,*)
         Write (*, '("Error(euler): matrix improper or not unitary")')
         Write (*, '(" Determinant : ", G18.10)') det
         Write (*,*)
         Stop
      End If
      If ((Abs(rot(3, 1)) .Gt. eps) .Or. (Abs(rot(3, 2)) .Gt. eps)) &
     & Then
         ang (1) = Atan2 (rot(3, 2), rot(3, 1))
         If (ang(1) .Lt. 0.d0) ang (1) = ang (1) + twopi
         If (Abs(rot(3, 1)) .Gt. Abs(rot(3, 2))) Then
            ang (2) = Atan2 (rot(3, 1)/Cos(ang(1)), rot(3, 3))
         Else
            ang (2) = Atan2 (rot(3, 2)/Sin(ang(1)), rot(3, 3))
         End If
         ang (3) = Atan2 (rot(2, 3),-rot(1, 3))
         If (ang(3) .Lt. 0.d0) ang (3) = ang (3) + twopi
      Else
         ang (1) = Atan2 (rot(1, 2), rot(1, 1))
         If (ang(1) .Lt. 0.d0) ang (1) = ang (1) + twopi
         If (rot(3, 3) .Gt. 0.d0) Then
            ang (2) = 0.d0
            ang (3) = 0.d0
         Else
            ang (2) = pi
            ang (3) = pi
         End If
      End If
      Return
End Subroutine
!EOC
