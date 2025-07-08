module rttddft_Wavefunction_test
  use constants, only: zone, zi, sqrt_two
  use exciting_mpi, only: mpiinfo
  use math_utils, only: all_close
  use precision, only: dp, i32
  use rttddft_Wavefunction, only: obtain_occupations
  use unit_test_framework, only : unit_test_type

  implicit none

  private

  public :: rttddft_Wavefunction_test_driver

  real(dp), parameter :: tol = 1.0e-10_dp
 
contains

  subroutine rttddft_Wavefunction_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional, intent(in) :: kill_on_failure
    
    type(unit_test_type) :: test_report
    integer(i32), parameter :: n_assertions_test_obtain_occupations = 2
    integer(i32), parameter :: n_assertions = n_assertions_test_obtain_occupations

    character(len=*), parameter :: module_tested = 'rttddft_Wavefunction'

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_obtain_occupations( 'SE', test_report )

    ! report results
    if ( present( kill_on_failure ) ) then
      call test_report%report( module_tested, kill_on_failure )
    else
      call test_report%report( module_tested )
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine

  
  subroutine test_obtain_occupations( method, test_report )
    !> Name of the propagator to be tested
    character(len=*), intent(in) :: method
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_identifier = "test_obtain_occupations"
    integer(i32), parameter :: n_states_gnd = 4
    integer(i32), parameter :: n_states = 3
    integer(i32), parameter :: n_kpt = 2
    real(dp), parameter :: occ_gnd(n_states_gnd, n_kpt) = reshape( &
                                    [2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
                                     2.0_dp, 1.8_dp, 0.2_dp, 0.0_dp], [n_states_gnd, n_kpt])
    real(dp), parameter :: occ_expected(n_states_gnd, n_kpt) = reshape( &
                                     [  1.5_dp,   1.5_dp, 0.5_dp, 0.5_dp, &
                                      1.928_dp, 1.872_dp, 0.1_dp, 0.1_dp], [n_states_gnd, n_kpt])
    complex(dp), parameter :: proj(n_states_gnd, n_states, n_kpt) = reshape(  &
                                     [    zone/sqrt_two,     zone/sqrt_two, (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), & !1st kpt
                                       (0.0_dp, 0.5_dp), -(0.0_dp, 0.5_dp), (0.5_dp, 0.0_dp), (0.5_dp, 0.0_dp), & !1st kpt
                                       (0.0_dp, 0.0_dp),  (0.0_dp, 0.0_dp), (0.0_dp, 0.8_dp), (0.0_dp, 0.6_dp), & !1st kpt
                                       (0.8_dp, 0.0_dp),  (0.6_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), & !2nd kpt
                                      -(0.0_dp, 0.6_dp),  (0.0_dp, 0.8_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), & !2nd kpt
                                       (0.0_dp, 0.0_dp),  (0.0_dp, 0.0_dp),      zi/sqrt_two,    zone/sqrt_two],& !2nd kpt
                                      [n_states_gnd, n_states, n_kpt] )
    real(dp), allocatable :: occ(:, :)                                                                  

    call obtain_occupations( proj, occ_gnd, occ )
    call test_report%assert( all( shape(occ) == shape(occ_expected) ), &
      message=test_identifier//' - wrong shape.')
    call test_report%assert( all_close( occ , occ_expected, tol=tol ), &
      message=test_identifier//' - calculated occ does not match reference.')
  end subroutine


end module