!> Unit tests for [[rttddft_pmat]]
module rttddft_pmat_test
  use constants, only: real_one, real_zero
  use exciting_mpi, only: mpiinfo
  use file_utils, only: delete_file
  use math_utils, only: all_close
  use mock_arrays, only: complex_hermitian_matrix_5x5, complex_matrix_5x5
  use mod_atoms, only: natmtot
  use modmpi, only: barrier
  use precision, only: dp, i32
  use rttddft_file_formats, only: file_handler, hdf5
  use rttddft_io_unformatted, only: delete_pmat_binary_file, delete_pmat_mt_binary_file
  use rttddft_pmat, only: pmat_set
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none

  private
  
  public :: rttddft_pmat_test_driver

  real(dp), parameter :: tol = 1.0e-7_dp

contains

  subroutine rttddft_pmat_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes if an assertion fails
    logical, optional, intent(in) :: kill_on_failure
    
    type(unit_test_type) :: test_report
    character(len=*), parameter :: module_tested = 'rttddft_pmat'

    ! Initialize test object
    call test_report%init( mpiglobal )

    ! Run and assert tests
    call test_pmat_read_write( mpiglobal, test_report )

    ! Report results
    call test_report%report( module_tested, kill_on_failure )

    ! Finalise test object
    call test_report%finalise()
  end subroutine

  !> Test the write and read subroutines of [[rttddft_pmat::pmat_set]] class
  subroutine test_pmat_read_write( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    integer(i32) :: p_dimension, first_kpt, last_kpt, test_counter
    integer(i32), parameter :: n_kpt_per_proc = 2
    character(len=*), parameter :: test_id = "test_pmat_read_write"
    character(len=:), allocatable :: test_case
    type(file_handler) :: hdf5_handler

    p_dimension = 5
    first_kpt = (mpiglobal%rank)*n_kpt_per_proc + 1
    last_kpt = first_kpt + n_kpt_per_proc - 1

    test_counter = 1
    test_case = 'hermitian matrix'
    call procedures_for_unit_test( complex_hermitian_matrix_5x5, is_hermitian=.true., use_pmat_MT=.false. )

    test_case = 'generic matrix'
    call procedures_for_unit_test( complex_matrix_5x5, is_hermitian=.false., use_pmat_MT=.false. )

    test_case = 'generic matrix and pmatMT'
    call procedures_for_unit_test( complex_matrix_5x5, is_hermitian=.false., use_pmat_MT=.true. )
#ifdef _HDF5_ 
    hdf5_handler%file_format = hdf5
    hdf5_handler%file_name = "rt.h5"
    hdf5_handler%path = "./"

    test_case = 'hermitian matrix HDF5'
    call procedures_for_unit_test( complex_hermitian_matrix_5x5, is_hermitian=.true., use_pmat_MT=.false., handler=hdf5_handler )

    test_case = 'generic matrix HDF5'
    call procedures_for_unit_test( complex_matrix_5x5, is_hermitian=.false., use_pmat_MT=.false., handler=hdf5_handler )
    
    test_case = 'generic matrix and pmatMT HDF5'
    call procedures_for_unit_test( complex_matrix_5x5, is_hermitian=.false., use_pmat_MT=.true., handler=hdf5_handler )
#endif
    contains
      !> (private) execute a series of common procedures for a unit test
      subroutine procedures_for_unit_test( matrix, is_hermitian, use_pmat_MT, handler )
        !> Matrix to be used in the unit test
        complex(dp), intent(in) :: matrix(:, :)
        !> Is the matrix hermitian?
        logical, intent(in) :: is_hermitian
        !> Use pmatMT?
        logical, intent(in) :: use_pmat_MT
        !> File handler
        type(file_handler), intent(in), optional :: handler

        integer(i32) :: i, j
        integer(i32), allocatable :: n_kpt
        integer(i32), parameter :: n_cartesian = 3
        logical :: use_hdf5
        complex(dp), allocatable :: px_ref(:, :, :), py_ref(:, :, :), pz_ref(:, :, :)
        complex(dp), allocatable :: pmat_MT_ref(:, :, :, :, :)
        type(pmat_set) :: pmat

        if( use_pmat_MT ) natmtot = 1
        call pmat%allocate( p_dimension, first_kpt, last_kpt, is_hermitian, use_pmat_MT, is_LAPWLO_basis=.true. )
        do i = first_kpt, last_kpt
          pmat%components(1)%array(:, :, i) = reshape( matrix*(1._dp + i*0.532_dp), [p_dimension, p_dimension] )
          pmat%components(2)%array(:, :, i) = reshape( matrix*(1._dp - i*0.532_dp), [p_dimension, p_dimension] )
          pmat%components(3)%array(:, :, i) = reshape( matrix, [p_dimension, p_dimension] )
        end do
        px_ref = pmat%components(1)%array
        py_ref = pmat%components(2)%array
        pz_ref = pmat%components(3)%array
        if( use_pmat_MT ) then
          do j = 1, natmtot
            pmat%MT(:, :, 1, j, :) = 0.73_dp*pmat%components(1)%array
            pmat%MT(:, :, 2, j, :) = 0.86_dp*pmat%components(2)%array
            pmat%MT(:, :, 3, j, :) = 0.19_dp*pmat%components(3)%array
          end do
          pmat_MT_ref = pmat%MT
        end if
        if( present(handler) ) allocate( n_kpt, source=(mpiglobal%procs*n_kpt_per_proc) )
        call pmat%write_to_file( mpiglobal, handler, n_kpt )
        call pmat%read_from_file( mpiglobal, handler )
        call test_report%assert( all_close( pmat%components(1)%array, px_ref, tol ), &
          wrapper_message( ' x' ) )
        call test_report%assert( all_close( pmat%components(2)%array, py_ref, tol ), &
          wrapper_message( ' y' ) )
        call test_report%assert( all_close( pmat%components(3)%array, pz_ref, tol ), &
          wrapper_message( ' z' ) )
        if( use_pmat_MT ) then
          do j = 1, n_cartesian
            call test_report%assert( all_close( pmat%MT(:, :, j, :, :), pmat_MT_ref(:, :, j, :, :), tol ), &
              wrapper_message( ' MT direction ' // to_char( j ) ) )
          end do
        end if
        call barrier()
        if( mpiglobal%is_root ) then
          use_hdf5 = .false.
          if( present(handler) ) use_hdf5 = ( handler%file_format == hdf5 )
          if( .not. use_hdf5 ) then
            call delete_pmat_binary_file( )
            if( use_pmat_MT ) call delete_pmat_MT_binary_file( )
          else
            call delete_file( handler%file_name, i )
          end if
        end if
        test_counter = test_counter + 1
      end subroutine

      !> (private) report message, given a test case and test number
      function wrapper_message( direction ) result(resulting_message)
        !> Direction: x, y or z
        character(len=*), intent(in) :: direction
        character(len=:), allocatable :: resulting_message
        resulting_message = report_message( test_id, test_case // direction, test_counter )
      end function
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