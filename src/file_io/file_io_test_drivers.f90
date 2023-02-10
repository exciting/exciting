module file_io_test_drivers
    use modmpi, only: mpiinfo
    ! Load test drivers here
    use os_utils_test, only: run_os_utils_test_driver
    use block_data_file_test, only: run_block_data_file_test_driver
    
    private
    public :: file_io_test_driver
  
  contains
    subroutine file_io_test_driver(mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(inout) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails 
      logical, optional :: kill_on_failure 
  
      ! Call test drivers here
      call run_os_utils_test_driver(mpiglobal, kill_on_failure)
      call run_block_data_file_test_driver(mpiglobal, kill_on_failure)
      
    end subroutine file_io_test_driver
  
end module file_io_test_drivers
