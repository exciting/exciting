!> Module with the unit tests related to groundstate
module groundstate_test_drivers
  use cdft_tests, only: run_cdft_test_driver
  use ghost_band_filter_test, only: run_ghost_band_filter_test_driver
  use general_find_vbm_cbm_test, only: general_find_vbm_cbm_test_driver
  use modmpi, only: mpiinfo

  implicit none
  
  private

  public :: groundstate_test_driver

contains

!> Call all unit tests related to groundstate
subroutine groundstate_test_driver(mpiglobal, kill_on_failure)
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal
  !> Kill the program before the test driver finishes
  !> if an assertion fails 
  logical, optional :: kill_on_failure 

  logical :: kill_if_fails

  kill_if_fails = .false.
  if( present(kill_on_failure) ) kill_if_fails = kill_on_failure
  
  ! Call test drivers here
  call run_cdft_test_driver( mpiglobal, kill_if_fails )
  call run_ghost_band_filter_test_driver(mpiglobal, kill_if_fails)
  call general_find_vbm_cbm_test_driver( mpiglobal, kill_if_fails )
end subroutine

end module
