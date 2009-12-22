!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sdelta_sq
! !INTERFACE:
Real (8) Function sdelta_sq (x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the square-wave pulse approximation to the Dirac delta function
!   $$ \tilde\delta(x)=\left\{\begin{array}{ll}
!    1 & \quad |x|\le 1/2 \\
!    0 & \quad |x|>1/2 \end{array}\right. $$
!
! !REVISION HISTORY:
!   Created July 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: x
      If (Abs(x) .Le. 0.5d0) Then
         sdelta_sq = 1.d0
      Else
         sdelta_sq = 0.d0
      End If
      Return
End Function
!EOC
