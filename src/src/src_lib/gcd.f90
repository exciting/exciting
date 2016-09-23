!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: gcd
! !INTERFACE:
Integer Function gcd (x, y)
! !INPUT/OUTPUT PARAMETERS:
!   x : first integer (in,integer)
!   y : second integer (in,integer)
! !DESCRIPTION:
!   Computes the greatest common divisor (GCD) of two integers using Euclid's
!   algorithm.
!
! !REVISION HISTORY:
!   Created September 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: x
      Integer, Intent (In) :: y
! local variables
      Integer :: a, b, c
      If (x .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(gcd): x <= 0 : ", I8)') x
         Write (*,*)
         Stop
      End If
      If (y .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(gcd): y <= 0 : ", I8)') y
         Write (*,*)
         Stop
      End If
      If (x .Ge. y) Then
         a = x
         b = y
      Else
         a = y
         b = x
      End If
10    Continue
      c = Mod (a, b)
      a = b
      b = c
      If (c .Gt. 0) Go To 10
      gcd = a
      Return
End Function
!EOC
