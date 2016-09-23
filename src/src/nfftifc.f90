!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: nfftifc
! !INTERFACE:
!
!
Subroutine nfftifc (n)
! !INPUT/OUTPUT PARAMETERS:
!   n : required/avalable grid size (in,integer)
! !DESCRIPTION:
!   Interface to the grid requirements of the fast Fourier transform routine.
!   Most routines restrict $n$ to specific prime factorisations. This routine
!   returns the next largest grid size allowed by the FFT routine.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (Inout) :: n
! local variables
      Integer :: i, j
! currently we use primes 2, 3 and 5
      Integer, Parameter :: np = 3
      Integer :: p (np)
      Data p / 2, 3, 5 /
      If (n .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(nfftifc): n <= 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
10    Continue
      i = n
      Do j = 1, np
         Do while (Mod(i, p(j)) .Eq.  0)
            i = i / p (j)
         End Do
      End Do
      If (i .Ne. 1) Then
         n = n + 1
         Go To 10
      End If
End Subroutine
!EOC
