!> Module for operations of the form \(\alpha x + y\)
module add_scaled_array
  use asserts, only: assert
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
  end interface

contains

  subroutine scaled_add_rank1_arrays_complex_dp(alpha, x, y)
    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: x(:)
    complex(dp), intent(inout) :: y(:)

    call assert( size(x) == size(y), "x and y must have same size" )
    call zaxpy( size(x), alpha, x, default_inc, y, default_inc )
  end subroutine

  subroutine scaled_add_rank2_arrays_complex_dp(alpha, x, y)
    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: x(:, :)
    complex(dp), intent(inout) :: y(:, :)

    call assert( all( shape(x) == shape(y) ), "x and y must have same shape" )
    call zaxpy( size(x), alpha, x, default_inc, y, default_inc )
  end subroutine

  subroutine scaled_add_rank3_arrays_complex_dp(alpha, x, y)
    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: x(:, :, :)
    complex(dp), intent(inout) :: y(:, :, :)

    call assert( all( shape(x) == shape(y) ), "x and y must have same shape" )
    call zaxpy( size(x), alpha, x, default_inc, y, default_inc )
  end subroutine
end module