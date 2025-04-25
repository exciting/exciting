!> Contains procedures to safely remap the pointer considering IFORT bug
module mod_pointer_remapping

    use precision, only: i32, dp

    implicit none

    private
    public :: remap_fortran_pointer

    interface remap_fortran_pointer
        module procedure remap_complex_dp_rank2
        module procedure remap_complex_dp_rank3
        module procedure remap_integer_i32_rank2
        module procedure remap_integer_i32_rank3
    end interface remap_fortran_pointer

contains

    ! Remap complex(dp) rank-2
    subroutine remap_complex_dp_rank2(ptr, lower_bounds, upper_bounds)
        complex(dp), pointer, contiguous, intent(inout) :: ptr(:,:)
        integer(i32), intent(in) :: lower_bounds(2)
        integer(i32), intent(in) :: upper_bounds(2)

#if defined(REMAPPING_BUG)
        complex(dp), pointer, contiguous :: temp_ptr(:,:)
        integer(i32) :: ptr_shape(2)
        ptr_shape = shape(ptr)
        temp_ptr(1:ptr_shape(1),1:ptr_shape(2)) => ptr(:,:)
        nullify(ptr)
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2)) => temp_ptr(:,:)
        nullify(temp_ptr)
#else
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2)) => ptr(:,:)
#endif
    end subroutine remap_complex_dp_rank2

    ! Remap complex(dp) rank-3
    subroutine remap_complex_dp_rank3(ptr, lower_bounds, upper_bounds)
        complex(dp), pointer, contiguous, intent(inout) :: ptr(:,:,:)
        integer(i32), intent(in) :: lower_bounds(3)
        integer(i32), intent(in) :: upper_bounds(3)

#if defined(REMAPPING_BUG)
        complex(dp), pointer, contiguous :: temp_ptr(:,:,:)
        integer(i32) :: ptr_shape(3)
        ptr_shape = shape(ptr)
        temp_ptr(1:ptr_shape(1),1:ptr_shape(2),1:ptr_shape(3)) => ptr(:,:,:)
        nullify(ptr)
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2), lower_bounds(3):upper_bounds(3)) => temp_ptr(:,:,:)
        nullify(temp_ptr)
#else
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2), lower_bounds(3):upper_bounds(3)) => ptr(:,:,:)
#endif
    end subroutine remap_complex_dp_rank3

    ! Remap integer(i32) rank-2
    subroutine remap_integer_i32_rank2(ptr, lower_bounds, upper_bounds)
        integer(i32), pointer, contiguous, intent(inout) :: ptr(:,:)
        integer(i32), intent(in) :: lower_bounds(2)
        integer(i32), intent(in) :: upper_bounds(2)

#if defined(REMAPPING_BUG)
        integer(i32), pointer, contiguous :: temp_ptr(:,:)
        integer(i32) :: ptr_shape(2)
        ptr_shape = shape(ptr)
        temp_ptr(1:ptr_shape(1),1:ptr_shape(2)) => ptr(:,:)
        nullify(ptr)
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2)) => temp_ptr(:,:)
        nullify(temp_ptr)
#else
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2)) => ptr(:,:)
#endif
    end subroutine remap_integer_i32_rank2

    ! Remap integer(i32) rank-3
    subroutine remap_integer_i32_rank3(ptr, lower_bounds, upper_bounds)
        integer(i32), pointer, contiguous, intent(inout) :: ptr(:,:,:)
        integer(i32), intent(in) :: lower_bounds(3)
        integer(i32), intent(in) :: upper_bounds(3)

#if defined(REMAPPING_BUG)
        integer(i32), pointer, contiguous :: temp_ptr(:,:,:)
        integer(i32) :: ptr_shape(3)
        ptr_shape = shape(ptr)
        temp_ptr(1:ptr_shape(1),1:ptr_shape(2),1:ptr_shape(3)) => ptr(:,:,:)
        nullify(ptr)
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2), lower_bounds(3):upper_bounds(3)) => temp_ptr(:,:,:)
        nullify(temp_ptr)
#else
        ptr(lower_bounds(1):upper_bounds(1), lower_bounds(2):upper_bounds(2), lower_bounds(3):upper_bounds(3)) => ptr(:,:,:)
#endif
    end subroutine remap_integer_i32_rank3

end module mod_pointer_remapping
