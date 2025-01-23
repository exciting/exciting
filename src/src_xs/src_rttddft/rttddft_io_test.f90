module rttddft_io_test
  use exciting_mpi, only: mpiinfo
  use math_utils, only: all_close
  use mock_arrays, only: real_matrix_5x7
  use modmpi, only: barrier
  use precision, only: dp, i32
  use rttddft_io, only: open_files_jpa, close_files_jpa, read_jpa, write_jpa, delete_jpa_files
  use rttddft_CurrentDensity, only: Current_Density_Field
  use rttddft_Polarization, only: Polarization
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none

  private

  public :: rttddft_io_test_driver

  real(dp), parameter :: tol = 1.0e-10_dp
 
contains

  subroutine rttddft_io_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional, intent(in) :: kill_on_failure
    
    type(unit_test_type) :: test_report
    integer(i32), parameter :: n_assertions_test_read_write_jpa = 18
    integer(i32), parameter :: n_assertions = n_assertions_test_read_write_jpa

    character(len=*), parameter :: module_tested = 'rttddft_io'

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_read_write_jpa( mpiglobal, test_report )

    ! report results
    if ( present( kill_on_failure ) ) then
      call test_report%report( module_tested, kill_on_failure )
    else
      call test_report%report( module_tested )
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine

  
  subroutine test_read_write_jpa( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_id = "test_read_write_jpa"
    character(len=*), parameter :: last = "last ", penult = "penultimate "
    character(len=*), parameter :: str_t = "time", str_j = "j_ind"
    character(len=*), parameter :: str_p = "p_vec", str_ai = "a_ind", str_at = "a_tot"
    character(len=*), parameter :: last_t = last // str_t, last_j = last // str_j
    character(len=*), parameter :: last_p = last // str_p, last_ai = last // str_ai, last_at = last // str_at
    character(len=*), parameter :: penult_j = penult // str_j, penult_p = penult // str_p
    character(len=*), parameter :: penult_ai = penult // str_ai, penult_at = penult // str_at
    real(dp) :: time
    real(dp), allocatable :: times(:), aux(:, :), aux2(:, :)
    type(Current_Density_Field) :: j_ind_t, j_ind_t_minus_dt
    type(Current_Density_Field), allocatable :: j_ind_list(:)
    type(Polarization) :: p_vec_t, p_vec_t_minus_dt
    type(Polarization), allocatable :: p_vec_list(:)
    type(Vector_Potential_Field) :: a_ind_t, a_tot_t, a_ind_t_minus_dt, a_tot_t_minus_dt
    type(Vector_Potential_Field), allocatable :: a_ind_list(:), a_tot_list(:)
    integer(i32) :: n, i, test_counter
    logical :: my_rank_writes

    my_rank_writes = mpiglobal%is_root
    n = 10
    times = [(i*1._dp, i = 1, n)]
    aux = reshape( real_matrix_5x7, [3, n] )
    aux2 = reshape( transpose(real_matrix_5x7), [3, n] )
    allocate( j_ind_list(n), p_vec_list(n), a_ind_list(n), a_tot_list(n) )
    do i = 1, n
      j_ind_list(i)%components = aux(:, i)
      p_vec_list(i)%components = aux2(:, i)
      a_ind_list(i)%components = aux(:, i) - aux2(:, i)**2
      a_tot_list(i)%components = aux2(:, i) + aux(:, i)**2
    end do
    if( my_rank_writes ) then
      call open_files_jpa( new=.true. )
      call write_jpa( times, j_ind_list )
      call write_jpa( times, p_vec_list )
      call write_jpa( times, a_ind_list, a_tot_list )
      call close_files_jpa
    end if
    call barrier

    test_counter = 1
    call read_jpa( time, j_ind_t )
    call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter ) )
    associate( j => j_ind_t%components, j_ref => j_ind_list(n)%components )
      call test_report%assert( all_close( j, j_ref, tol ), report_message( test_id, last_j, test_counter) )
      
      test_counter = 2
      call read_jpa( time, j_ind_t, j_ind_t_minus_dt )
      call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
      call test_report%assert( all_close( j, j_ref, tol ), report_message( test_id, last_j, test_counter) )
    end associate
    associate( j => j_ind_t_minus_dt%components, j_ref => j_ind_list(n-1)%components )
      call test_report%assert( all_close( j, j_ref, tol ), report_message( test_id, penult_j, test_counter) )
    end associate
    
    test_counter = 3
    call read_jpa( time, p_vec_t )
    call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
    associate( p => p_vec_t%components, p_ref => p_vec_list(n)%components )
      call test_report%assert( all_close( p, p_ref, tol ), report_message( test_id, last_p, test_counter) )

      test_counter = 4
      call read_jpa( time, p_vec_t, p_vec_t_minus_dt )
      call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
      call test_report%assert( all_close( p, p_ref, tol ), report_message( test_id, last_p, test_counter) )
    end associate
    associate( p => p_vec_t_minus_dt%components, p_ref => p_vec_list(n-1)%components )
      call test_report%assert( all_close( p, p_ref, tol ), report_message( test_id, penult_p, test_counter) )
    end associate

    test_counter = 5
    call read_jpa( time, a_ind_t, a_tot_t=a_tot_t )
    call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
    associate( ai => a_ind_t%components, ai_ref => a_ind_list(n)%components, at => a_tot_t%components, at_ref => a_tot_list(n)%components )
      call test_report%assert( all_close( ai , ai_ref, tol ), report_message( test_id, last_ai, test_counter) )
      call test_report%assert( all_close( at , at_ref, tol ), report_message( test_id, last_at, test_counter) )

      test_counter = 6
      call read_jpa( time, a_ind_t, a_ind_t_minus_dt, a_tot_t, a_tot_t_minus_dt )
      call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
      call test_report%assert( all_close( ai , ai_ref, tol ), report_message( test_id, last_ai, test_counter) )
      call test_report%assert( all_close( at , at_ref, tol ), report_message( test_id, last_at, test_counter) )
    end associate
    associate( ai => a_ind_t_minus_dt%components, ai_ref => a_ind_list(n-1)%components, at => a_tot_t_minus_dt%components, at_ref => a_tot_list(n-1)%components )
      call test_report%assert( all_close( ai , ai_ref, tol ), report_message( test_id, penult_ai, test_counter) )
      call test_report%assert( all_close( at , at_ref, tol ), report_message( test_id, penult_at, test_counter) )
    end associate

    call barrier
    if( my_rank_writes ) call delete_jpa_files
  end subroutine

  !> (private) generates a report message, given a test case and test number
  function report_message( test_id, test_case, test_number ) result(message)
    !> Name of the subroutine calling `report_message`
    character(len=*), intent(in) :: test_id
    !> Test (sub)identifier, pointing out which case in [[test_id]] is tested
    character(len=*), intent(in) :: test_case
    !> Number for the test
    integer(i32), intent(in) :: test_number
    !> Message to return
    character(len=:), allocatable :: message
    message = test_id // ' - ' // test_case // ' does not match reference: test - ' // to_char(test_number)
  end function


end module