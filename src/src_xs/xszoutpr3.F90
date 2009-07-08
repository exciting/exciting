

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

module m_xszoutpr3
  implicit none
contains

!BOP
! !ROUTINE: xszoutpr3
! !INTERFACE:


subroutine xszoutpr3(n1, n2, alpha, x, y, a)
! !INPUT/OUTPUT PARAMETERS:
!   n1,n2 : size of vectors and matrix, respectively (in,integer)
!   alpha : complex constant (in,complex)
!   x     : first input vector (in,complex(n1))
!   y     : second input vector (in,complex(n2))
!   a     : output matrix (out,complex(n1,n2))
! !DESCRIPTION:
!   Performs the rank-2 operation
!   $$ A_{ij}\rightarrow\alpha{\bf x}_i{\bf y}_j+A_{ij}. $$
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: n1, n2
    complex(8), intent(in) :: alpha, x(:), y(:)
    complex(8), intent(inout) :: a(:, :)
    ! local variables
    integer :: i2
    complex(8) :: yt
    do i2=1, n2
       yt=y(i2)
       call zaxpy(n1, alpha*yt, x, 1, a(1, i2), 1)
    end do
  end subroutine xszoutpr3
!EOC

end module m_xszoutpr3
