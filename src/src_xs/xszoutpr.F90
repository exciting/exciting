
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

module m_xszoutpr
  implicit none
contains

  subroutine xszoutpr(n1,n2,alpha,x,y,a)
    !
    ! perform the operation: A_ij -> A_ij + alpha*conjg(x_i)*y_j
    !
    implicit none
    ! arguments
    integer, intent(in) :: n1,n2
    complex(8), intent(in) :: alpha, x(:), y(:)
    complex(8), intent(inout) :: a(:,:)
    ! local variables
    integer :: i2
    ! allocatable arrays
    complex(8), allocatable :: h(:)
    allocate(h(size(x)))
    h=conjg(x)
    do i2=1,n2
       call zaxpy(n1,alpha*y(i2),h,1,a(1,i2),1)
    end do
    deallocate(h)
  end subroutine xszoutpr

end module m_xszoutpr
