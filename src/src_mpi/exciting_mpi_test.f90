!> Module for unit tests for for the functions in exciting_mpi.
module exciting_mpi_test
  use exciting_mpi, only: xmpi_allgather, xmpi_allgatherv
  use math_utils, only: all_close
  use modmpi, only: mpiinfo, distribute_loop
  use mock_arrays, only: fill_array, value_map_complex_rank1, &
    value_map_complex_rank2, value_map_complex_rank3
  use precision, only: dp, i32, sp
  use unit_test_framework, only: unit_test_type

  implicit none

  private
  public :: exciting_mpi_test_driver
  character(len=*), parameter :: module_tested = "exciting_mpi"
  real(dp), parameter :: tol = 1.e-14_dp

contains

  !> Run tests for modmpi
  subroutine exciting_mpi_test_driver( mpiglobal, kill_on_failure )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report

    call test_report%init( mpiglobal )

    ! Run unit tests
    call test_xmpi_allgather( test_report, mpiglobal)

    if ( present( kill_on_failure ) ) then
      call test_report%report( module_tested, kill_on_failure )
    else
      call test_report%report( module_tested )
    end if

    call test_report%finalise()
  end subroutine exciting_mpi_test_driver

  subroutine test_xmpi_allgather( test_report, mpiglobal )
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal

    character(len=*), parameter :: test_identifier = "test_xmpi_allgather"
    integer(i32), parameter :: n_elements = 17, n_rows = 31, n_columns = 7
    integer(i32) :: i, first, last, n_local
    complex(dp) :: complex_dp_rank_1(n_elements), complex_dp_rank_1_ref(n_elements), &
      complex_dp_rank_2(n_rows, n_elements), complex_dp_rank_2_ref(n_rows, n_elements), &
      complex_dp_rank_3(n_rows, n_columns, n_elements), complex_dp_rank_3_ref(n_rows, n_columns, n_elements)
    real(dp) :: real_dp_rank_1(n_elements), real_dp_rank_1_ref(n_elements), &
      real_dp_rank_2(n_rows, n_elements), real_dp_rank_2_ref(n_rows, n_elements), &
      real_dp_rank_3(n_rows, n_columns, n_elements), real_dp_rank_3_ref(n_rows, n_columns, n_elements)
    real(sp) :: real_sp_rank_1(n_elements), real_sp_rank_1_ref(n_elements), &
      real_sp_rank_2(n_rows, n_elements), real_sp_rank_2_ref(n_rows, n_elements), &
      real_sp_rank_3(n_rows, n_columns, n_elements), real_sp_rank_3_ref(n_rows, n_columns, n_elements)
    integer(i32) :: integer_i32_rank_1(n_elements), integer_i32_rank_1_ref(n_elements), &
      integer_i32_rank_2(n_rows, n_elements), integer_i32_rank_2_ref(n_rows, n_elements), &
      integer_i32_rank_3(n_rows, n_columns, n_elements), integer_i32_rank_3_ref(n_rows, n_columns, n_elements)
    integer(i32), allocatable :: receive_buffer(:), receive_buffer_ref(:)

    call distribute_loop( mpiglobal, n_elements, first, last )
    n_local = last - first + 1

    call fill_array( complex_dp_rank_1_ref, value_map_complex_rank1 )
    call fill_array( complex_dp_rank_2_ref, value_map_complex_rank2 )
    call fill_array( complex_dp_rank_3_ref, value_map_complex_rank3 )

    complex_dp_rank_1(first : last) = complex_dp_rank_1_ref(first : last)
    complex_dp_rank_2(:, first : last) = complex_dp_rank_2_ref(:, first : last)
    complex_dp_rank_3(:, :, first : last) = complex_dp_rank_3_ref(:, :, first : last)

    real_dp_rank_1 = real( complex_dp_rank_1, kind = dp )
    real_dp_rank_1_ref = real( complex_dp_rank_1_ref, kind = dp )
    real_dp_rank_2 = real( complex_dp_rank_2, kind = dp )
    real_dp_rank_2_ref = real( complex_dp_rank_2_ref, kind = dp )
    real_dp_rank_3 = real( complex_dp_rank_3, kind = dp )
    real_dp_rank_3_ref = real( complex_dp_rank_3_ref, kind = dp )

    real_sp_rank_1 = real( complex_dp_rank_1, kind = sp )
    real_sp_rank_1_ref = real( complex_dp_rank_1_ref, kind = sp )
    real_sp_rank_2 = real( complex_dp_rank_2, kind = sp )
    real_sp_rank_2_ref = real( complex_dp_rank_2_ref, kind = sp )
    real_sp_rank_3 = real( complex_dp_rank_3, kind = sp )
    real_sp_rank_3_ref = real( complex_dp_rank_3_ref, kind = sp )

    integer_i32_rank_1 = int( real_dp_rank_1, kind = i32 )
    integer_i32_rank_1_ref = real( real_dp_rank_1_ref, kind = i32 )
    integer_i32_rank_2 = real( real_dp_rank_2, kind = i32 )
    integer_i32_rank_2_ref = real( real_dp_rank_2_ref, kind = i32 )
    integer_i32_rank_3 = real( real_dp_rank_3, kind = i32 )
    integer_i32_rank_3_ref = real( real_dp_rank_3_ref, kind = i32 )

    call xmpi_allgatherv( mpiglobal, complex_dp_rank_1, n_local )
    call test_report%assert( all_close( complex_dp_rank_1, complex_dp_rank_1_ref, tol = tol ), &
      test_identifier//"_complex_dp_rank_1" )
    call xmpi_allgatherv( mpiglobal, complex_dp_rank_2, size(complex_dp_rank_2, 1) * n_local )
    call test_report%assert( all_close( complex_dp_rank_2, complex_dp_rank_2_ref, tol = tol ), &
      test_identifier//"_complex_dp_rank_2" )
    call xmpi_allgatherv( mpiglobal, complex_dp_rank_3, size(complex_dp_rank_3, 1) * size(complex_dp_rank_3, 2) * n_local )
    call test_report%assert( all_close( complex_dp_rank_3, complex_dp_rank_3_ref, tol = tol ), &
      test_identifier//"_complex_dp_rank_3" )

    call xmpi_allgatherv( mpiglobal, real_dp_rank_1, n_local )
    call test_report%assert( all_close( real_dp_rank_1, real_dp_rank_1_ref, tol = tol ), &
      test_identifier//"_real_dp_rank_1" )
    call xmpi_allgatherv( mpiglobal, real_dp_rank_2, size(real_dp_rank_2, 1) * n_local )
    call test_report%assert( all_close( real_dp_rank_2, real_dp_rank_2_ref, tol = tol ), &
      test_identifier//"_real_dp_rank_2" )
    call xmpi_allgatherv( mpiglobal, real_dp_rank_3, size(real_dp_rank_3, 1) * size(real_dp_rank_3, 2) * n_local )
    call test_report%assert( all_close( real_dp_rank_3, real_dp_rank_3_ref, tol = tol ), &
      test_identifier//"_real_dp_rank_3" )

    call xmpi_allgatherv( mpiglobal, real_sp_rank_1, n_local )
    call test_report%assert( all_close( real_sp_rank_1, real_sp_rank_1_ref, tol = real( tol, kind = sp ) ), &
      test_identifier//"_real_sp_rank_1" )
    call xmpi_allgatherv( mpiglobal, real_sp_rank_2, size(real_sp_rank_2, 1) * n_local )
    call test_report%assert( all_close( real_sp_rank_2, real_sp_rank_2_ref, tol = real( tol, kind = sp ) ), &
      test_identifier//"_real_sp_rank_2" )
    call xmpi_allgatherv( mpiglobal, real_sp_rank_3, size(real_sp_rank_3, 1) * size(real_sp_rank_3, 2) * n_local )
    call test_report%assert( all_close( real_sp_rank_3, real_sp_rank_3_ref, tol = real( tol, kind = sp ) ), &
      test_identifier//"_real_sp_rank_3" )

    call xmpi_allgatherv( mpiglobal, integer_i32_rank_1, n_local )
    call test_report%assert( all( integer_i32_rank_1 == integer_i32_rank_1_ref ), &
      test_identifier//"_integer_i32_rank_1" )
    call xmpi_allgatherv( mpiglobal, integer_i32_rank_2, size(integer_i32_rank_2, 1) * n_local )
    call test_report%assert( all( integer_i32_rank_2 == integer_i32_rank_2_ref ), &
      test_identifier//"_integer_i32_rank_2" )
    call xmpi_allgatherv( mpiglobal, integer_i32_rank_3, size(integer_i32_rank_3, 1) * size(integer_i32_rank_3, 2) * n_local )
    call test_report%assert( all(  integer_i32_rank_3 == integer_i32_rank_3_ref ), &
      test_identifier//"_integer_i32_rank_3" )

    allocate( receive_buffer_ref(mpiglobal%procs) )
    do i = 1, size( receive_buffer_ref )
      receive_buffer_ref(i) = i
    end do
    call xmpi_allgather( mpiglobal, mpiglobal%rank + 1, receive_buffer )
    call test_report%assert( all( receive_buffer == receive_buffer_ref ), &
      test_identifier//"_integer_i32_not_in_place" )
  end subroutine test_xmpi_allgather
   
end module exciting_mpi_test

