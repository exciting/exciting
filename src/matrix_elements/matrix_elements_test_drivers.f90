module matrix_elements_test_drivers
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use matrix_elements_test, &
    only: me_test_init, &
          me_test_orthogonality, &
          me_test_mt, &
          me_test_sf, &
          me_test_grad, &
          me_test_sfgrad, &
          me_test_operator, &
          me_test_gauss

  implicit none
  private

  public :: matrix_elements_test_driver

  contains

    subroutine matrix_elements_test_driver(mpiglobal, kill_on_failure)
      use matrix_elements_test, only: ng
      !> mpi information
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails 
      logical, optional :: kill_on_failure 
  
      !> Test object
      type(unit_test_type) :: test_report
      
      !> Number of assertions
      integer, parameter :: n_assertions = 1+1+3+3*(3+3)+3*3*2+ng+3*ng

      call test_report%init( n_assertions, mpiglobal)

      ! initialize tests
      call me_test_init
      ! perform orthogonality test
      ! 1 assertion
      call me_test_orthogonality( test_report)
      ! perform test of muffin-tin volume integrals
      ! 1 assertion
      call me_test_mt( test_report)
      ! perform tests of muffin-tin surface integrals
      ! 3 assertions
      call me_test_sf( test_report)
      ! perform tests of muffin-tin volume integrals with gradients
      ! wave functions
      ! 18 assertions
      call me_test_grad( test_report)
      ! perform tests of muffin-tin surface integrals with gradients
      ! wave functions
      ! 18 assertions
      call me_test_sfgrad( test_report)
      ! perform orthogonality tests for different plane-wave operators
      ! ng assertions
      call me_test_operator( test_report)
      ! perform tests on Gauss' theorem
      ! 3*ng assertions
      call me_test_gauss( test_report)
      
      if( present(kill_on_failure)) then
        call test_report%report( 'matrix_elements', kill_on_failure)
      else
        call test_report%report( 'matrix_elements')
      end if
   
      call test_report%finalise()
    end subroutine matrix_elements_test_driver

end module matrix_elements_test_drivers
