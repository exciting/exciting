!> Module for collecting unit test drivers for modules in the fft directory.
module fft_test_drivers
  use modmpi, only: mpiinfo
  ! Load test drivers here
  !use fftw_mpi_wrapper_test, only: fftw_mpi_wrapper_test_driver
  private
  public :: fft_test_driver

  contains

  subroutine fft_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure 

    ! Call test drivers here
    !call fftw_mpi_wrapper_test_driver(mpiglobal, kill_on_failure)
  end subroutine fft_test_driver

end module fft_test_drivers