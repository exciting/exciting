!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: stheta_sq
! !INTERFACE:
Real (8) Function stheta_sq (x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the Heaviside step function corresponding to the square-wave pulse
!   approximation to the Dirac delta function
!   $$ \tilde\Theta(x)=\left\{\begin{array}{ll}
!    0 & \quad x \le -1/2 \\
!    x+1/2 & \quad -1/2 < x < 1/2 \\
!    1 & \quad x\ge 1 \end{array}\right. $$
!
! !REVISION HISTORY:
!   Created July 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: x
      If (x .Le.-0.5d0) Then
         stheta_sq = 0.d0
         Return
      End If
      If (x .Lt. 0.5d0) Then
         stheta_sq = x + 0.5d0
      Else
         stheta_sq = 1.d0
      End If
      Return
End Function
!EOC
