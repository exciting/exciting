!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sdelta_fd
! !INTERFACE:
Real (8) Function sdelta_fd (x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the Fermi-Dirac approximation to the Dirac delta function
!   $$ \tilde\delta(x)=\frac{e^{-x}}{(1+e^{-x})^2}. $$
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: x
! local variables
      Real (8) :: t1
      If (Abs(x) .Gt. 50.d0) Then
         sdelta_fd = 0.d0
         Return
      End If
      t1 = Exp (-x)
      sdelta_fd = t1 / ((1.d0+t1)**2)
      Return
End Function
!EOC
