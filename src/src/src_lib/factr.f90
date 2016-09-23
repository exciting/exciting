!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: factr
! !INTERFACE:
Real (8) Function factr (n, d)
! !INPUT/OUTPUT PARAMETERS:
!   n : numerator (in,integer)
!   d : denominator (in,integer)
! !DESCRIPTION:
!   Returns the ratio $n!/d!$ for $n,d\ge 0$. Performs no under- or overflow
!   checking.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Integer, Intent (In) :: d
! local variables
      Integer :: i
      If (n .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(factr): n < 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (d .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(factr): d < 0 : ", I8)') d
         Write (*,*)
         Stop
      End If
      If (n .Lt. d) Then
         factr = 1.d0 / dble (n+1)
         Do i = n + 2, d
            factr = factr / dble (i)
         End Do
      Else If (n .Eq. d) Then
         factr = 1.d0
      Else
         factr = dble (d+1)
         Do i = d + 2, n
            factr = factr * dble (i)
         End Do
      End If
      Return
End Function
!EOC
