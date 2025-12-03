! Created by  on 24/03/2022.
!> Test [[to_char_conversion]]
module string_utils_test
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type

  use string_utils, only: split, validate_filename

  implicit none

  private
  public :: string_utils_test_driver

  contains

   !> Run tests for the lattice module
  subroutine string_utils_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> test object
    type(unit_test_type) :: test_report

    call test_report%init( mpiglobal)

    call test_split(test_report)

    call test_validate_filename(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('string_utils', kill_on_failure)
    else
      call test_report%report('string_utils')
    end if
  end subroutine string_utils_test_driver

  !> Test `[[split]]`.
  subroutine test_split(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report

    character(*), parameter :: string = 'ADOijfafd.FASPOIPWR.SARWQsd-abc-asdegesd'
    character(:), allocatable :: first, second, delimiter 

    delimiter = '.'
    call split(string, delimiter, first, second)

    call test_report%assert(first == 'ADOijfafd', &
            'test_split: expected all characters before the first period as first.')

    call test_report%assert(second == 'FASPOIPWR.SARWQsd-abc-asdegesd', &
            'test_split: expected all characters after the first period as first.')

  end subroutine test_split

  !> Test `[[validate_filename]]`.
  subroutine test_validate_filename(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report

    character(:), allocatable :: filename

    ! Valid file name
    filename = 'asijhfdadf_sadiom-asfgs.h5'
    call test_report%assert(validate_filename(filename, '.h5'), &
            'test_validate_filename: Expected true for valid filename.')
    call test_report%assert(validate_filename(filename, 'h5'), &
            'test_validate_filename: Expected true for valid filename.')

    ! Non-valid file name (two periods before the ending)
    filename = 'asijhfdadf_sadiom..h5'
    call test_report%assert(.not. validate_filename(filename, '.h5'), &
            'test_validate_filename: Expected false for non-valid filename.')
    call test_report%assert(.not. validate_filename(filename, 'h5'), &
            'test_validate_filename: Expected false for non-valid filename.')

    ! Non-valid file name (forbidden character in filename body)
    filename = 'asijhfdadf%sadiom.h5'
    call test_report%assert(.not. validate_filename(filename, '.h5'), &
            'test_validate_filename: Expected false for non-valid filename.')
    call test_report%assert(.not. validate_filename(filename, 'h5'), &
            'test_validate_filename: Expected false for non-valid filename.')

    ! Non-valid file name (forbidden character in filename ending)
    filename = 'asijhfdadf_sadiom.h^5'
    call test_report%assert(.not. validate_filename(filename, '.h^5'), &
            'test_validate_filename: Expected false for non-valid filename.')
    call test_report%assert(.not. validate_filename(filename, 'h^5'), &
            'test_validate_filename: Expected false for non-valid filename.')

    ! Non-valid file name (forbidden character in filename body & ending)
    filename = 'asijhfdadf#sadiom.h)5'
    call test_report%assert(.not. validate_filename(filename, '.h)5'), &
            'test_validate_filename: Expected false for non-valid filename.')
    call test_report%assert(.not. validate_filename(filename, 'h)5'), &
            'test_validate_filename: Expected false for non-valid filename.')                   

  end subroutine test_validate_filename

end module string_utils_test