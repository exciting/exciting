!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: stheta_fd
! !INTERFACE:
Real (8) Function stheta_fd (x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the Fermi-Dirac approximation to the Heaviside step function
!   $$ \tilde\Theta(x)=\frac{1}{1+e^{-x}}. $$
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: x
      If (x .Gt. 50.d0) Then
         stheta_fd = 1.d0
         Return
      End If
      If (x .Lt.-50.d0) Then
         stheta_fd = 0.d0
         Return
      End If
      stheta_fd = 1.d0 / (1.d0+Exp(-x))
      Return
End Function
!EOC
