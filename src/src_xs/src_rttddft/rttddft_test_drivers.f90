module rttddft_test_drivers
  use modmpi, only: mpiinfo
  use propagators_test, only: propagators_test_driver
  use rttddft_io_test, only: rttddft_io_test_driver
  use rttddft_Wavefunction_test, only: rttddft_Wavefunction_test_driver

  implicit none
  
  private

  public :: rttddft_test_driver

  logical, parameter :: kill_on_failure_default = .true.

contains
  subroutine rttddft_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional, intent(in) :: kill_on_failure 

    logical :: kill_on_failure_

    kill_on_failure_ = kill_on_failure_default
    if( present(kill_on_failure) ) kill_on_failure_ = kill_on_failure

    ! Call test drivers here
    call propagators_test_driver( mpiglobal, kill_on_failure_ )
    call rttddft_io_test_driver( mpiglobal, kill_on_failure_ )
    call rttddft_Wavefunction_test_driver( mpiglobal, kill_on_failure_ )
  end subroutine
end module
