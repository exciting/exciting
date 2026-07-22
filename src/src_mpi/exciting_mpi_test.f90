!> Module for unit tests for for the functions in exciting_mpi.
module exciting_mpi_test
  use exciting_mpi, only: xmpi_allgather, xmpi_allgatherv, xmpi_bcast, xmpi_reduce, xmpi_allreduce, xmpi_gather
  use math_utils, only: all_close
  use modmpi, only: mpiinfo, distribute_loop
  use mock_arrays, only: fill_array, value_map_complex_rank1, &
    value_map_complex_rank2, value_map_complex_rank3
  use precision, only: dp, i32, long_int, sp
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
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report

    call test_report%init( mpiglobal )

    ! Run unit tests
    call test_xmpi_allgather( test_report, mpiglobal )
    call test_xmpi_allgatherv_edge_cases( test_report, mpiglobal )
    call test_xmpi_bcast( test_report, mpiglobal )
    call test_xmpi_reduce( test_report, mpiglobal )
    call test_xmpi_allreduce( test_report, mpiglobal )
    call test_xmpi_gather( test_report, mpiglobal )

    if ( present( kill_on_failure ) ) then
      call test_report%report( module_tested, kill_on_failure )
    else
      call test_report%report( module_tested )
    end if

    call test_report%finalise()
  end subroutine exciting_mpi_test_driver

subroutine test_xmpi_bcast( test_report, mpiglobal )
    type(unit_test_type), intent(inout) :: test_report
    type(mpiinfo), intent(in) :: mpiglobal

    character(len=*), parameter :: test_identifier = "test_xmpi_bcast"
    integer(i32), parameter :: n_elements = 10, n_rows = 5
    integer(i32) :: alt_root
    
    real(dp) :: real_dp_rank_1(n_elements), real_dp_rank_1_ref(n_elements)
    integer(i32) :: int_i32_rank_2(n_rows, n_elements), int_i32_rank_2_ref(n_rows, n_elements)

    real_dp_rank_1_ref = 3.14159_dp
    int_i32_rank_2_ref = 42

    if ( mpiglobal%is_root ) then
      real_dp_rank_1 = real_dp_rank_1_ref
      int_i32_rank_2 = int_i32_rank_2_ref
    else
      real_dp_rank_1 = 0.0_dp
      int_i32_rank_2 = 0
    end if

    call xmpi_bcast( mpiglobal, real_dp_rank_1 )
    call xmpi_bcast( mpiglobal, int_i32_rank_2 )

    call test_report%assert( all_close( real_dp_rank_1, real_dp_rank_1_ref, tol = tol ), &
      test_identifier//"_real_dp_rank_1" )
    call test_report%assert( all( int_i32_rank_2 == int_i32_rank_2_ref ), &
      test_identifier//"_integer_i32_rank_2" )

    alt_root = mpiglobal%procs - 1

    if ( mpiglobal%rank == alt_root ) then
      real_dp_rank_1 = real_dp_rank_1_ref * 2.0_dp
      int_i32_rank_2 = int_i32_rank_2_ref * 2
    else
      real_dp_rank_1 = 0.0_dp
      int_i32_rank_2 = 0
    end if

    call xmpi_bcast( mpiglobal, real_dp_rank_1, bcasting_rank=alt_root )
    call xmpi_bcast( mpiglobal, int_i32_rank_2, bcasting_rank=alt_root )

    call test_report%assert( all_close( real_dp_rank_1, real_dp_rank_1_ref * 2.0_dp, tol = tol ), &
      test_identifier//"_alt_root_real_dp_rank_1" )
    call test_report%assert( all( int_i32_rank_2 == int_i32_rank_2_ref * 2 ), &
      test_identifier//"_alt_root_integer_i32_rank_2" )

  end subroutine test_xmpi_bcast

  subroutine test_xmpi_reduce( test_report, mpiglobal )
    type(unit_test_type), intent(inout) :: test_report
    type(mpiinfo), intent(in) :: mpiglobal

    character(len=*), parameter :: test_identifier = "test_xmpi_reduce"
    integer(i32), parameter :: n_elements = 15
    real(dp) :: real_dp_rank_1(n_elements), real_dp_rank_1_ref(n_elements)

    real_dp_rank_1 = 1.0_dp
    real_dp_rank_1_ref = 1.0_dp * mpiglobal%procs

    call xmpi_reduce( real_dp_rank_1, mpiglobal )

    if ( mpiglobal%is_root ) then
      call test_report%assert( all_close( real_dp_rank_1, real_dp_rank_1_ref, tol = tol ), &
        test_identifier//"_real_dp_rank_1" )
    end if
  end subroutine test_xmpi_reduce

  subroutine test_xmpi_allreduce( test_report, mpiglobal )
    type(unit_test_type), intent(inout) :: test_report
    type(mpiinfo), intent(in) :: mpiglobal

    character(len=*), parameter :: test_identifier = "test_xmpi_allreduce"
    integer(i32), parameter :: n_elements = 15
    complex(dp) :: cmplx_dp_rank_1(n_elements), cmplx_dp_rank_1_ref(n_elements)

    cmplx_dp_rank_1 = cmplx(1.0_dp, 1.0_dp, kind=dp)
    cmplx_dp_rank_1_ref = cmplx(1.0_dp * mpiglobal%procs, 1.0_dp * mpiglobal%procs, kind=dp)

    call xmpi_allreduce( cmplx_dp_rank_1, mpiglobal )

    call test_report%assert( all_close( cmplx_dp_rank_1, cmplx_dp_rank_1_ref, tol = tol ), &
      test_identifier//"_complex_dp_rank_1" )
  end subroutine test_xmpi_allreduce

  subroutine test_xmpi_gather( test_report, mpiglobal )
    type(unit_test_type), intent(inout) :: test_report
    type(mpiinfo), intent(in) :: mpiglobal

    character(len=*), parameter :: test_identifier = "test_xmpi_gather"
    integer(i32), parameter :: n_local = 5
    integer(i32) :: i
    integer(i32) :: send_arr(n_local)
    integer(i32), allocatable :: recv_arr(:), recv_arr_ref(:)

    send_arr = mpiglobal%rank + 1

    call xmpi_gather( mpiglobal, send_arr, recv_arr )

    if ( mpiglobal%is_root ) then
      allocate( recv_arr_ref(n_local * mpiglobal%procs) )
      do i = 1, mpiglobal%procs
        recv_arr_ref((i-1)*n_local + 1 : i*n_local) = i
      end do

      call test_report%assert( allocated(recv_arr), test_identifier//"_is_allocated" )
      call test_report%assert( size(recv_arr) == size(recv_arr_ref), test_identifier//"_size_check" )
      call test_report%assert( all(recv_arr == recv_arr_ref), test_identifier//"_data_check" )
    end if
  end subroutine test_xmpi_gather

  subroutine test_xmpi_allgather( test_report, mpiglobal )
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI environment
    type(mpiinfo), intent(in) :: mpiglobal

    character(len=*), parameter :: test_identifier = "test_xmpi_allgather"
    integer(i32), parameter :: n_elements = 17, n_rows = 31, n_columns = 7
    integer(i32) :: i, first, last, n_local
    integer(long_int) :: n_local_large
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
    integer_i32_rank_1_ref = int( real_dp_rank_1_ref, kind = i32 )
    integer_i32_rank_2 = int( real_dp_rank_2, kind = i32 )
    integer_i32_rank_2_ref = int( real_dp_rank_2_ref, kind = i32 )
    integer_i32_rank_3 = int( real_dp_rank_3, kind = i32 )
    integer_i32_rank_3_ref = int( real_dp_rank_3_ref, kind = i32 )

    n_local_large = int( n_local, kind = long_int )
    call xmpi_allgatherv( mpiglobal, complex_dp_rank_1, n_local_large )
    call test_report%assert( all_close( complex_dp_rank_1, complex_dp_rank_1_ref, tol = tol ), &
      test_identifier//"_complex_dp_rank_1_large" )

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

  subroutine test_xmpi_allgatherv_edge_cases( test_report, mpiglobal )
    type(unit_test_type), intent(inout) :: test_report
    type(mpiinfo), intent(in) :: mpiglobal

    character(len=*), parameter :: test_identifier = "test_xmpi_allgatherv_edges"
    integer(i32) :: chunk_size, total_elements, i, displ_start
    integer(long_int) :: chunk_size_large
    real(dp), allocatable :: buffer(:), buffer_ref(:)

    ! --- TEST 1: Non-uniform distribution ---
    ! Rank `r` contributes `r + 1` elements.
    ! This heavily stresses the displacement logic and the `flat_buf` slicing workaround.
    chunk_size = mpiglobal%rank + 1
    total_elements = mpiglobal%procs * (mpiglobal%procs + 1) / 2
    chunk_size_large = int(chunk_size, kind=long_int)

    allocate(buffer(total_elements), buffer_ref(total_elements))
    buffer = 0.0_dp

    ! Build reference array (1, 2, 3... N)
    do i = 1, total_elements
      buffer_ref(i) = real(i, dp)
    end do

    ! Simulate MPI_IN_PLACE: Place local data exactly where it belongs globally
    displ_start = mpiglobal%rank * (mpiglobal%rank + 1) / 2
    buffer(displ_start + 1 : displ_start + chunk_size) = buffer_ref(displ_start + 1 : displ_start + chunk_size)

    call xmpi_allgatherv( mpiglobal, buffer, chunk_size_large )
    call test_report%assert( all_close( buffer, buffer_ref, tol = tol ), &
      test_identifier//"_non_uniform" )

    ! --- TEST 2: Zero-contribution edge case ---
    ! Rank 0 has all the data, everyone else has 0. 
    ! This ensures the workaround's `max(1, chunk_size)` safeguard prevents crashes.
    buffer = 0.0_dp
    if ( mpiglobal%rank == 0 ) then
      chunk_size = total_elements
      buffer = buffer_ref
    else
      chunk_size = 0
    end if
    chunk_size_large = int(chunk_size, kind=long_int)

    call xmpi_allgatherv( mpiglobal, buffer, chunk_size_large )
    call test_report%assert( all_close( buffer, buffer_ref, tol = tol ), &
      test_identifier//"_zero_chunk" )

    deallocate(buffer, buffer_ref)
  end subroutine test_xmpi_allgatherv_edge_cases
   
end module exciting_mpi_test