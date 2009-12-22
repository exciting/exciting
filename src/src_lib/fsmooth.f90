!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: fsmooth
!
!
Subroutine fsmooth (m, n, ld, f)
! !INPUT/OUTPUT PARAMETERS:
!   m  : number of 3-point running averages to perform (in,integer)
!   n  : number of point (in,integer)
!   ld : leading dimension (in,integer)
!   f  : function array (inout,real(ld,n))
! !DESCRIPTION:
!   Removes numerical noise from a function by performing $m$ successive
!   3-point running averages on the data. The endpoints are kept fixed.
!
! !REVISION HISTORY:
!   Created December 2005 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: m
      Integer, Intent (In) :: n
      Integer, Intent (In) :: ld
      Real (8), Intent (Inout) :: f (ld, n)
! local variables
      Integer :: i, j
! automatic arrays
      Real (8) :: g (n)
      Do i = 1, m
         Do j = 2, n - 1
            g (j) = 0.3333333333333333333d0 * (f(1, j-1)+f(1, j)+f(1, &
           & j+1))
         End Do
         f (1, 2:n-1) = g (2:n-1)
      End Do
      Return
End Subroutine
!EOC
