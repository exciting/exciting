module rttddft_arrays_utils
  use multi_index_conversion, only: indices_to_composite_index
  use precision, only: i32, dp

  implicit none
  
  private

  public :: map_array_to_pointer

  contains
  !> Remap indexes of a rank-two array to a rank-one pointer as illustrated in
  !> `ptr(new_first:new_last) => arr(1:size_dim1, first_dim2:last_dim2)`
  subroutine remap_indexes( size_dim1, first_dim2, last_dim2, new_first, new_last )
    integer(i32), intent(in) :: size_dim1, first_dim2, last_dim2
    integer(i32), intent(out) :: new_first, new_last

    new_first = indices_to_composite_index( [1, first_dim2], [size_dim1, last_dim2] )
    new_last = indices_to_composite_index( [size_dim1, last_dim2], [size_dim1, last_dim2] )
  end subroutine

  !> Map a rank-5 array to rank-4 pointer
  subroutine map_array_to_pointer( first, array, ptr_rank4 )
    integer(i32), intent(in) :: first
    complex(dp), contiguous, target, intent(in) :: array(:, :, :, :, first:)
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)
    ! Local variables
    integer(i32) :: a, b, c, d, last, first_rk4, last_rk4

    a = size( array, 1 )
    b = size( array, 2 )
    c = size( array, 3 )
    d = size( array, 4 )
    last = ubound( array, 5 )
    call remap_indexes( d, first, last, first_rk4, last_rk4)
    ptr_rank4(1:a,1:b,1:c,first_rk4:last_rk4) => array
  end subroutine

end module