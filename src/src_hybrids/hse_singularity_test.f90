!> Module for unit tests for the functions in the hse_singularity module.

module hse_singularity_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use math_utils, only: all_close
  use unit_test_framework, only : unit_test_type
  use hse_singularity , only : hse_singularity_exact_solution, & 
                             & hse_singularity_Taylor_expansion    
  implicit none
  private
  public :: hse_singularity_test_driver

contains

  !> Run tests for funtions to compute the hse_singularity
  subroutine hse_singularity_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> Test report object
    type(unit_test_type) :: test_report

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Run and assert tests

    call test_hse_singularity_exact(test_report)

    call test_hse_singularity_Taylor(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('hse_singularity', kill_on_failure)
    else
      call test_report%report('hse_singularity')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine hse_singularity_test_driver 


  !> Test all_close for real and complex arrays
  subroutine test_hse_singularity_exact(test_report)

    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    !> Real test input
    real(dp) :: omega_hse 
    real(dp) :: volume_unit_cell
    
    !> Integer test input
    integer :: number_k_points

    !> Real test reference 
    real(dp) :: hse_integral_reference

    !> Real test calculated value
    real(dp) :: hse_integral
    
    !> Tolerance of deviation for tests       
    real(dp), parameter :: tol = 1.e-10_dp    
   
    volume_unit_cell = 848.21_dp
    omega_hse = 0.11_dp  
    number_k_points = 8
    hse_integral_reference = 203.53123301172664_dp
    
    hse_integral= hse_singularity_exact_solution(omega_hse,volume_unit_cell, number_k_points)
    call test_report%assert(all_close(hse_integral, hse_integral_reference,tol), &
                  & 'Expected solution for the hse integral &
                  & computed exactly is 203.53123301172664_dp')

  end subroutine test_hse_singularity_exact

  subroutine test_hse_singularity_Taylor(test_report)

    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    !> Real test input
    real(dp) :: omega_hse 
    
    !> Real test reference 
    real(dp) :: hse_integral_reference

    !> Real test calculated value
    real(dp) :: hse_integral
    
    !> Tolerance of deviation for tests    
    real(dp), parameter :: tol = 1.e-10_dp 
    
    omega_hse = 0.11_dp
    hse_integral_reference = 259.6357564950242_dp
    
    hse_integral= hse_singularity_Taylor_expansion(omega_hse)
    call test_report%assert(all_close(hse_integral,hse_integral_reference,tol), &
                  & 'Expected solution for the hse integral computed with & 
                  & the Taylor expandsion method is 259.6357564950242_dp')

  end subroutine test_hse_singularity_Taylor


end module

