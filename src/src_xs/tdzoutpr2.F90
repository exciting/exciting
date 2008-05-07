
module m_tdzoutpr2
  implicit none
contains

  subroutine tdzoutpr2(n1,n2,alpha,x,y,a)
    !
    ! perform the operation: A_ij -> A_ij + alpha*x_i*conjg(y_j)
    !
    implicit none
    ! arguments
    integer, intent(in) :: n1,n2
    complex(8), intent(in) :: alpha, x(:), y(:)
    complex(8), intent(inout) :: a(:,:)
    ! local variables
    integer :: i2
    complex(8) :: yc,h(size(x))

    h=x
    do i2 = 1, n2
       yc=conjg(y(i2))
       call zaxpy(n1,alpha*yc,h,1,a(1,i2),1)
    end do

  end subroutine tdzoutpr2

end module m_tdzoutpr2

module m_tdzoutpr3
  implicit none
contains

  subroutine tdzoutpr3(n1,n2,alpha,x,y,a)
    !
    ! perform the operation: A_ij -> A_ij + alpha*x_i*y_j
    !
    implicit none
    ! arguments
    integer, intent(in) :: n1,n2
    complex(8), intent(in) :: alpha, x(:), y(:)
    complex(8), intent(inout) :: a(:,:)
    ! local variables
    integer :: i2
    complex(8) :: yc,h(size(x))

    h=x
    do i2 = 1, n2
       yc=y(i2)
       call zaxpy(n1,alpha*yc,h,1,a(1,i2),1)
    end do

  end subroutine tdzoutpr3

end module m_tdzoutpr3
