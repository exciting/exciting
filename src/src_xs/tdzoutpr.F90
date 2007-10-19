
module m_tdzoutpr
  implicit none
contains

  subroutine tdzoutpr(n1,n2,alpha,x,y,a)
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
    ! automatic arrays
    complex(8) :: h(size(x))

    h = conjg(x)
    do i2 = 1, n2
       call zaxpy(n1,alpha*y(i2),h,1,a(1,i2),1)
    end do
    
  end subroutine tdzoutpr

end module m_tdzoutpr
