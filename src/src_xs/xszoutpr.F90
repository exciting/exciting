!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
Module m_xszoutpr
      Implicit None
Contains
!
!BOP
! !ROUTINE: xszoutpr
! !INTERFACE:
!
!
      Subroutine xszoutpr (n1, n2, alpha, x, y, a)
! !INPUT/OUTPUT PARAMETERS:
!   n1,n2 : size of vectors and matrix, respectively (in,integer)
!   alpha : complex constant (in,complex)
!   x     : first input vector (in,complex(n1))
!   y     : second input vector (in,complex(n2))
!   a     : output matrix (out,complex(n1,n2))
! !DESCRIPTION:
!   Performs the rank-2 operation
!   $$ A_{ij}\rightarrow\alpha{\bf x}_i^*{\bf y}_j+A_{ij}. $$
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: n1, n2
         Complex (8), Intent (In) :: alpha, x (:), y (:)
         Complex (8), Intent (Inout) :: a (:, :)
    ! local variables
         Integer :: i2
    ! allocatable arrays
         Complex (8), Allocatable :: h (:)
         Allocate (h(size(x)))
         h = conjg (x)
         Do i2 = 1, n2
            Call zaxpy (n1, alpha*y(i2), h, 1, a(1, i2), 1)
         End Do
         Deallocate (h)
      End Subroutine xszoutpr
!EOC
!
End Module m_xszoutpr
