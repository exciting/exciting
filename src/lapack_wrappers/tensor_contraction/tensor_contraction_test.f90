!> Unit tests for general matrix multiplication.
module tensor_contraction_test
    use precision, only: dp
    use constants, only: zi
    use modmpi, only: mpiinfo
    use unit_test_framework, only: unit_test_type
    use math_utils, only: all_close
    use mock_arrays, only: fill_array, &
                           value_map_real_rank2, &
                           value_map_real_rank3, &
                           value_map_complex_rank2, &
                           value_map_complex_rank3, &
                           value_map_complex_rank4
    use tensor_contractions, only: complex_tensor_contraction_dp,&
                                   real_tensor_contraction_dp,&
                                   shape_rank_2
    
    implicit none

    private
    public :: tensor_contraction_test_driver

contains

    !> Runs tests for tensor_contraction.
    subroutine tensor_contraction_test_driver(mpiglobal, kill_on_failure)
        !> mpi environment
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program upon failure of an assertion
        logical, intent(in), optional :: kill_on_failure

        !> Test report object
        type(unit_test_type) :: test_report
        !> Number of assertions
        integer, parameter :: n_assertions = 10

        call test_report%init(n_assertions, mpiglobal)

        ! Run unit tests
        call test_complex_tensor_contraction_rank3_rank2_to_rank3(test_report)

        call test_real_tensor_contraction_rank3_rank2_to_rank3(test_report)

        call test_shape_rank_2(test_report)


        if (present(kill_on_failure)) then
            call test_report%report('tensor_contraction', kill_on_failure)
        else
            call test_report%report('tensor_contraction')
        end if

        call test_report%finalise()
    end subroutine tensor_contraction_test_driver

    
    !> Test contraction over one index of real tensors of
    !> rank-3 and rank-2, respectively. Yields a complex rank-3 tensor.
    subroutine test_real_tensor_contraction_rank3_rank2_to_rank3(test_report)
        use constants, only: zzero
        !> Test object
        type(unit_test_type), intent(inout) :: test_report

        real(dp), allocatable :: A(:, :, :), B(:, :), &
                                    C(:, :, :), C_ref(:, :, :)
        integer :: i, j, k, l, dims_C(3), size_contracted_dim

        dims_C = (/5, 5, 5/)

        ! A and B not transposed
        allocate (A(5, 5, 3))
        allocate (B(3, 5))
        allocate (C(dims_C(1), dims_C(2), dims_C(3)))
        allocate (C_ref(5, 5, 5))

        size_contracted_dim = size(A, dim=3)

        C = 0._dp
        C_ref = 0._dp

        call fill_array(A, value_map_real_rank3)
        call fill_array(B, value_map_real_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(i, j, l)*B(l, k)
                    end do
                end do
            end do
        end do

        call real_tensor_contraction_dp(A, shape(A), B, shape(B), C, shape(C))

        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: A * B')

        ! A and B transposed
        deallocate (A, B)
        allocate (A(3, 5, 5))
        allocate (B(5, 3))

        size_contracted_dim = size(A, dim=1)

        C = 0._dp
        C_ref = 0._dp

        call fill_array(A, value_map_real_rank3)
        call fill_array(B, value_map_real_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(l, i, j)*B(k, l)
                    end do
                end do
            end do
        end do

        call real_tensor_contraction_dp(A, shape(A), B, shape(B),&
                                         C, shape(C), trans_A='t', trans_B='t')
        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: transpose(A) * transpose(B)')

        ! Only A transposed
        deallocate (A, B)
        allocate (A(3, 5, 5))
        allocate (B(3, 5))

        size_contracted_dim = size(A, dim=1)

        C = 0._dp
        C_ref = 0._dp

        call fill_array(A, value_map_real_rank3)
        call fill_array(B, value_map_real_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(l, i, j)*B(l, k)
                    end do
                end do
            end do
        end do

        call real_tensor_contraction_dp(A, shape(A), B, shape(B),&
                                        C, shape(C), trans_A='t', trans_B='n')
        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: transpose(A) * B')

        ! Only B transposed
        deallocate (A, B)
        allocate (A(5, 5, 3))
        allocate (B(5, 3))

        size_contracted_dim = size(A, dim=3)

        C = 0._dp
        C_ref = 0._dp

        call fill_array(A, value_map_real_rank3)
        call fill_array(B, value_map_real_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(i, j, l)*B(k, l)
                    end do
                end do
            end do
        end do

        call real_tensor_contraction_dp(A, shape(A), B, shape(B),&
                                         C, shape(C), trans_A='n', trans_B='t')
        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: A * transpose(B)')

    end subroutine test_real_tensor_contraction_rank3_rank2_to_rank3
    
    
    !> Test contraction over one index of complex tensors of
    !> rank-3 and rank-2, respectively. Yields a complex rank-3 tensor.
    subroutine test_complex_tensor_contraction_rank3_rank2_to_rank3(test_report)
        use constants, only: zzero
        !> Test object
        type(unit_test_type), intent(inout) :: test_report

        complex(dp), allocatable :: A(:, :, :), B(:, :), &
                                    C(:, :, :), C_ref(:, :, :)
        integer :: i, j, k, l, dims_C(3), size_contracted_dim

        dims_C = (/5, 5, 5/)

        ! A and B not transposed
        allocate (A(5, 5, 3))
        allocate (B(3, 5))
        allocate (C(dims_C(1), dims_C(2), dims_C(3)))
        allocate (C_ref(5, 5, 5))

        size_contracted_dim = size(A, dim=3)

        C = zzero
        C_ref = zzero

        call fill_array(A, value_map_complex_rank3)
        call fill_array(B, value_map_complex_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(i, j, l)*B(l, k)
                    end do
                end do
            end do
        end do

        call complex_tensor_contraction_dp(A, shape(A), B, shape(B), C, shape(C))

        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: A * B')

        ! A and B transposed
        deallocate (A, B)
        allocate (A(3, 5, 5))
        allocate (B(5, 3))

        size_contracted_dim = size(A, dim=1)

        C = zzero
        C_ref = zzero

        call fill_array(A, value_map_complex_rank3)
        call fill_array(B, value_map_complex_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(l, i, j)*B(k, l)
                    end do
                end do
            end do
        end do

        call complex_tensor_contraction_dp(A, shape(A), B, shape(B),&
                                         C, shape(C), trans_A='t', trans_B='t')
        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: transpose(A) * transpose(B)')

        ! Only A transposed
        deallocate (A, B)
        allocate (A(3, 5, 5))
        allocate (B(3, 5))

        size_contracted_dim = size(A, dim=1)

        C = zzero
        C_ref = zzero

        call fill_array(A, value_map_complex_rank3)
        call fill_array(B, value_map_complex_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(l, i, j)*B(l, k)
                    end do
                end do
            end do
        end do

        call complex_tensor_contraction_dp(A, shape(A), B, shape(B),&
                                        C, shape(C), trans_A='t', trans_B='n')
        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: transpose(A) * B')

        ! Only B transposed
        deallocate (A, B)
        allocate (A(5, 5, 3))
        allocate (B(5, 3))

        size_contracted_dim = size(A, dim=3)

        C = zzero
        C_ref = zzero

        call fill_array(A, value_map_complex_rank3)
        call fill_array(B, value_map_complex_rank2)

        do i = 1, dims_C(1)
            do j = 1, dims_C(2)
                do k = 1, dims_C(3)
                    do l = 1, size_contracted_dim
                        C_ref(i, j, k) = C_ref(i, j, k) + A(i, j, l)*B(k, l)
                    end do
                end do
            end do
        end do

        call complex_tensor_contraction_dp(A, shape(A), B, shape(B),&
                                         C, shape(C), trans_A='n', trans_B='t')
        call test_report%assert(all_close(C, C_ref), &
                                'Test tensor contraction over one 1 index for  complex tensors: Rank-3 * Rank-2 -> Rank-3 . &
                                Expected: A * transpose(B)')

    end subroutine test_complex_tensor_contraction_rank3_rank2_to_rank3

    !> For tensors with rank > 2, tests that the corresponding two-dimensional shape
    !> used for the tensor contractions is correct. 
    !> For tensor A with shape(A) = [1,2,3,4,5], it is tested that for a contraction
    !> over two dimensions the two-dimensional shape is given by 
    !> shape_A_2d = [2,3*4*5] if the contracted dimensions are the leading and by
    !> shape_A_2d = [1*2*3,4*5] if the contracted dimensions are the trailing ones.

    subroutine test_shape_rank_2(test_report)
        use constants, only: zzero
        !> Test object
        type(unit_test_type), intent(inout) :: test_report
        !> Test tensor
        complex(dp) :: A(1,2,3,4,5)
        !> Number of contracted dimensions
        integer, parameter :: n_contracted = 2
        !> Two-dimensional shape of A used for contraction
        integer :: shape_A_2d(2)


        shape_A_2d = shape_rank_2(shape(A), n_contracted, .false.)

        call test_report%assert(all(shape_A_2d==[1*2,3*4*5]), &
                                'Test shape_rank_2 for a rank-5 array and 2 &
                                contracted dimensions (leading dimensions).&
                                Expected: shape_A_2d=[1*2,3*4*5]')

        shape_A_2d = shape_rank_2(shape(A), n_contracted, .true.)

        call test_report%assert(all(shape_A_2d==[1*2*3,4*5]), &
                                'Test shape_rank_2 for a rank-5 array and 2 &
                                contracted dimensions (trailing dimensions).&
                                Expected: shape_A_2d=[1*2*3,4*5]')


    end subroutine

end module tensor_contraction_test
