
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zoutpr
! !INTERFACE:
subroutine zoutpr(n1,n2,alpha,x,y,a)
! !INPUT/OUTPUT PARAMETERS:
!   n1,n2 : size of vectors and matrix, respectively (in,integer)
!   alpha : complex constant (in,complex)
!   x     : first input vector (in,complex(n))
!   y     : second input vector (in,complex(n))
!   a     : output matrix
! !DESCRIPTION:
!   Performs the rank-2 operation
!   $$ A_{ij}\rightarrow\alpha{\bf x}_i^*{\bf y}_j+A_{ij}. $$
!
! !REVISION HISTORY:
!   Created April 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: n1,n2
  complex(8), intent(in) :: alpha, x(n1), y(n2)
  complex(8), intent(inout) :: a(n1,n2)
  ! local variables
  integer :: i2
  complex(8), allocatable :: h(:)
  allocate(h(size(x)))
  h(:)=conjg(x(:))
  do i2=1,n2
     call zaxpy(n1,alpha*y(i2),h,1,a(1,i2),1)
  end do
  deallocate(h)
end subroutine zoutpr
!EOC
