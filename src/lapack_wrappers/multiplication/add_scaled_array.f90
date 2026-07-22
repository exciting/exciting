!> Module for operations of the form \(\alpha x + y\)
module add_scaled_array
#include "asserts.fpp"
  use precision, only: dp, i32
  use lapack_f95_interfaces, only: zaxpy

  implicit none

  private

  public :: scaled_add

  integer(i32), parameter :: default_inc = 1

  !> Wrapper to `zaxpy`: \(y := \alpha x + y\)
  interface scaled_add
    module procedure scaled_add_rank1_arrays_complex_dp
    module procedure scaled_add_rank2_arrays_complex_dp
    module procedure scaled_add_rank3_arrays_complex_dp
    module procedure scaled_add_rank3_arrays_complex_dp_alpha_real_dp
  end interface

contains

  subroutine scaled_add_rank1_arrays_complex_dp(alpha, x, y)
    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: x(:)
    complex(dp), intent(inout) :: y(:)

    CALL_ASSERT( size(x) == size(y), "x and y must have same size" )
    call zaxpy( size(x), alpha, x, default_inc, y, default_inc )
  end subroutine

  subroutine scaled_add_rank2_arrays_complex_dp(alpha, x, y)
    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: x(:, :)
    complex(dp), intent(inout) :: y(:, :)

    CALL_ASSERT( all( shape(x) == shape(y) ), "x and y must have same shape" )
    call zaxpy( size(x), alpha, x, default_inc, y, default_inc )
  end subroutine

  subroutine scaled_add_rank3_arrays_complex_dp(alpha, x, y)
    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: x(:, :, :)
    complex(dp), intent(inout) :: y(:, :, :)

    CALL_ASSERT( all( shape(x) == shape(y) ), "x and y must have same shape" )
    call zaxpy( size(x), alpha, x, default_inc, y, default_inc )
  end subroutine

  subroutine scaled_add_rank3_arrays_complex_dp_alpha_real_dp(alpha, x, y)
    real(dp), intent(in) :: alpha
    complex(dp), intent(in) :: x(:, :, :)
    complex(dp), intent(inout) :: y(:, :, :)

    CALL_ASSERT( all( shape(x) == shape(y) ), "x and y must have same shape" )
    ! In this case, the best option is to use Fortran array operations
    y = y + alpha*x
  end subroutine
end module