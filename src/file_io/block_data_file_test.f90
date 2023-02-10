!> Module for unit tests for the functions in [[block_data_file(module)]].
module block_data_file_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use block_data_file

  implicit none
  private

  character(*), parameter :: TEST_DIR = 'BLOCK_DATA_FILE_TEST_DIR/'
  character(*), parameter :: TEST_FILE = 'BLOCK_DATA_FILE_TEST_FILE'
  integer, parameter :: data_block_dim(3) = [1, 2, 3]
  integer, parameter :: num_data_block = 3

  type(block_data_file_type) :: ifile, rfile, zfile

  public :: run_block_data_file_test_driver

  contains

    !> Run tests for [[block_data_file(module)]]
    subroutine run_block_data_file_test_driver(mpiglobal, kill_on_failure)
      use os_utils, only: make_directory, remove_directory
      !> mpi information
      type(mpiinfo), intent(inout) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails
      logical, optional :: kill_on_failure
      !> Test report object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 14 + 3*num_data_block

      ! Initialize test object
      call test_report%init(n_assertions, mpiglobal)

      ! create test directory (test covered by `os_utils_test`)
      call test_report%assert( make_directory(TEST_DIR, mpiglobal) == 0, &
        'Expected: Returned error code 0 from `make_directory`.')
      
      ! Run and assert tests
      call test_gen_block_data_file(test_report, mpiglobal)
      call test_write_and_read_block_data_file(test_report)
      call test_del_block_data_file(test_report, mpiglobal)

      ! remove test directory (test covered by `os_utils_test`)
      call test_report%assert( remove_directory(TEST_DIR, mpiglobal) == 0, &
        'Expected: Returned error code 0 from `remove_directory`.')

      ! report results
      if (present(kill_on_failure)) then
        call test_report%report('file_io_utils', kill_on_failure)
      else
        call test_report%report('file_io_utils')
      end if

      ! Finalise test object
      call test_report%finalise()
    end subroutine run_block_data_file_test_driver

    !> Test `block_data_file_gen` and `block_data_file_open`
    subroutine test_gen_block_data_file(test_report, comm)
      use os_utils, only: path_exists
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(inout) :: comm

      integer :: ierr

      ! create file objects and open files
      ifile = block_data_file_type( &
        trim(TEST_DIR)//trim(TEST_FILE)//'_INTEGER.BIN', data_block_dim, 0)
      call ifile%open(comm, delete_existing=.true.)
      call test_report%assert( path_exists(ifile%get_path(), ierr), &
        'Integer type file does not exists.')
      call test_report%assert( ifile%is_open(), &
        'Integer type file `ifile` is not open.')

      rfile = block_data_file_type( &
        trim(TEST_DIR)//trim(TEST_FILE)//'_DOUBLE.BIN', data_block_dim, 0._dp)
      call rfile%open(comm, delete_existing=.true.)
      call test_report%assert( path_exists(rfile%get_path(), ierr), &
        'Double type file does not exists.')
      call test_report%assert( rfile%is_open(), &
        'Double type file `rfile` is not open.')

      zfile = block_data_file_type( &
        trim(TEST_DIR)//trim(TEST_FILE)//'_DOUBLE_COMPLEX.BIN', data_block_dim, cmplx(0,0,dp))
      call zfile%open(comm, delete_existing=.true.)
      call test_report%assert( path_exists(zfile%get_path(), ierr), &
        'Double complex type file does not exists.')
      call test_report%assert( zfile%is_open(), &
        'Double complex type file `zfile` is not open.')
    end subroutine test_gen_block_data_file

    !> Test `block_data_file_close` and `block_data_file_delete`
    subroutine test_del_block_data_file(test_report, comm)
      use os_utils, only: path_exists
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(inout) :: comm

      integer :: ierr
      logical :: exists

      ! close and delete files
      call ifile%close(comm)
      call test_report%assert( .not. ifile%is_open(), &
        'Integer type file `ifile` is not closed.')
      call ifile%delete(comm)
      call test_report%assert( .not. path_exists(ifile%get_path(), ierr), &
        'Integer type file still exists.')

      call rfile%close(comm)
      call test_report%assert( .not. rfile%is_open(), &
        'Double type file `rfile` is not closed.')
      call rfile%delete(comm)
      call test_report%assert( .not. path_exists(rfile%get_path(), ierr), &
        'Double type file still exists.')

      call zfile%close(comm)
      call test_report%assert( .not. zfile%is_open(), &
        'Double complex type file is not closed.')
      call zfile%delete(comm)
      call test_report%assert( .not. path_exists(zfile%get_path(), ierr), &
        'Double complex type file still exists.')
    end subroutine test_del_block_data_file

    !> Test `block_data_file_read`
    subroutine test_write_and_read_block_data_file(test_report)
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: i, j, n

      integer, allocatable :: iblock_write(:,:,:,:), iblock_read(:,:,:), iblock1d(:)
      real(dp), allocatable :: rblock_write(:,:,:,:), rblock_read(:,:,:), rblock1d(:)
      complex(dp), allocatable :: zblock_write(:,:,:,:), zblock_read(:,:,:), zblock1d(:)

      n = product( data_block_dim )

      allocate( iblock1d(n), rblock1d(n), zblock1d(n) )
      allocate( iblock_write(data_block_dim(1), data_block_dim(2), data_block_dim(3), num_data_block) )
      allocate( rblock_write(data_block_dim(1), data_block_dim(2), data_block_dim(3), num_data_block) )
      allocate( zblock_write(data_block_dim(1), data_block_dim(2), data_block_dim(3), num_data_block) )
      allocate( iblock_read, source=reshape( iblock1d, data_block_dim ) )
      allocate( rblock_read, source=reshape( rblock1d, data_block_dim ) )
      allocate( zblock_read, source=reshape( zblock1d, data_block_dim ) )

      ! write data
      do i = 1, num_data_block
        do j = 1, n
          iblock1d(j) = i*j
          rblock1d(j) = 1._dp*i*j
          zblock1d(j) = cmplx( 1._dp*i*j, -1._dp*i*j, dp )
        end do

        iblock_write(:,:,:,i) = reshape( iblock1d, data_block_dim )
        call ifile%write( i, iblock_write(:,:,:,i) )

        rblock_write(:,:,:,i) = reshape( rblock1d, data_block_dim )
        call rfile%write( i, rblock_write(:,:,:,i) )

        zblock_write(:,:,:,i) = reshape( zblock1d, data_block_dim )
        call zfile%write( i, zblock_write(:,:,:,i) )
      end do

      ! read data
      do i = 1, num_data_block
        call ifile%read( i, iblock_read )
        call test_report%assert( all( iblock_read == iblock_write(:,:,:,i) ), &
          'Written and read values do not match for integer type file `ifile`.')

        call rfile%read( i, rblock_read )
        call test_report%assert( all( rblock_read == rblock_write(:,:,:,i) ), &
          'Written and read values do not match for double type file `rfile`.')

        call zfile%read( i, zblock_read )
        call test_report%assert( all( zblock_read == zblock_write(:,:,:,i) ), &
          'Written and read values do not match for double complex type file `zfile`.')
      end do
    end subroutine test_write_and_read_block_data_file

end module block_data_file_test
