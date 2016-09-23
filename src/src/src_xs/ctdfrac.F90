!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_ctdfrac
      Implicit None
Contains
!
!BOP
! !ROUTINE: ctdfrac
! !INTERFACE:
!
!
      Subroutine ctdfrac (n, a, b, f)
! !INPUT/OUTPUT PARAMETERS:
!   n     : depth of continued fraction (in,integer)
!   a     : a-coefficients (in,complex(n))
!   b     : b-coefficients (in,complex(0:n))
!   f     : continued fraction result
! !DESCRIPTION:
!   Straight forward evaluation of a continued fraction with depth $n$
!   $$ b_0+\cfrac{a_1}{b_1+\cfrac{a_2}
!   {b_2+\cfrac{a_3}{\cdots+\cfrac{a_n}{b_n}}}}  $$
!   without any checking of convergence.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: n
         Complex (8), Intent (In) :: a (n), b (0:n)
         Complex (8), Intent (Out) :: f
    ! local variables
         Integer :: j
         f = b (n)
         Do j = n, 1, - 1
            f = b (j-1) + a (j) / f
         End Do
      End Subroutine
!
End Module
!EOC
