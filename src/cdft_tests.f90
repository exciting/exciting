!> Module that contains the CDFT unit tests
module cdft_tests
  use constants, only: zi, zone, zzero
  use cdft, only: ExcitonCoefficients, occupy_cdft
  use math_utils, only: all_close
  use modmpi, only: mpiinfo
  use precision, only: i32, dp
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none

  private

  public :: run_cdft_test_driver

contains 
!> Run the CDFT unit tests
subroutine run_cdft_test_driver( mpiglobal, kill_on_failure )
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal
  !> Kill the program before the test driver finishes if an assertion fails
  logical, optional :: kill_on_failure
  
  type(unit_test_type) :: test_report
  integer, parameter :: n_assertions_test_get_ExcitonCoefficients_from_file = 5
  integer, parameter :: n_assertions_test_occupy_cdft = 2
  integer, parameter :: n_assertions = n_assertions_test_get_ExcitonCoefficients_from_file + &
                                       n_assertions_test_occupy_cdft 
  character(len=*), parameter :: test_driver_name = "cdft"

  call test_report%init( n_assertions, mpiglobal )

  ! Run and assert tests
  call test_ExcitonCoefficients_get_from_file( test_report, mpiglobal )
  call test_occupy_cdft( test_report )
  
  if (present(kill_on_failure)) then
    call test_report%report( test_driver_name, kill_on_failure )
  else
    call test_report%report( test_driver_name )
  end if

  call test_report%finalise()

end subroutine


!> Unit tests for [[cdft(module):ExcitonCoefficients_get_from_file(subroutine)]]
subroutine test_ExcitonCoefficients_get_from_file( test_report, mpiglobal )
  !> Unit test report
  type(unit_test_type) :: test_report
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal

  integer(i32), parameter :: n_non_zero_ref = 4
  integer(i32), parameter :: idx_vb_ref(n_non_zero_ref) = [4, 6, 8, 9] 
  integer(i32), parameter :: idx_cb_ref(n_non_zero_ref) = [11, 12, 11, 14]
  integer(i32), parameter :: idx_kpt_ref(n_non_zero_ref) = [1, 1, 7, 9] 
  character(len=*), parameter :: fake_name = "test"
  complex(dp), parameter :: weights_ref(n_non_zero_ref) = [ zzero, zi, zone, zi + zone ]
  real(dp), parameter :: tol = 1e-8_dp
  integer(i32) :: i, unit
  character(len=:), allocatable :: file_name
  integer(i32), allocatable :: idx_vb(:), idx_cb(:), idx_kpt(:)
  complex(dp), allocatable :: coeffs(:)
  type(ExcitonCoefficients) :: exc

  ! Each MPI rank reads/writes its own file
  file_name = fake_name // to_char( mpiglobal%rank )
  ! Write to file
  open( newunit = unit, file = trim(file_name), action = 'write' )
  write( unit, * ) n_non_zero_ref
  do i = 1, n_non_zero_ref
    write( unit, '(3I6, 2F18.10)' ) idx_vb_ref(i), idx_cb_ref(i), idx_kpt_ref(i), weights_ref(i)
  end do
  close( unit )

  call exc%get_from_file( file_name )
  call exc%get_attributes( idx_kpt, idx_vb, idx_cb, coeffs )
  call test_report%assert( n_non_zero_ref == size( coeffs ), "Wrong number of coefficients" )
  call test_report%assert( all( idx_kpt == idx_kpt_ref ), "k-point indexes do not match" )
  call test_report%assert( all( idx_vb == idx_vb_ref ), "valence-band indexes do not match" )
  call test_report%assert( all( idx_cb == idx_cb_ref ), "conduction-band indexes do not match" )
  call test_report%assert( all_close(coeffs, weights_ref, tol), "exciton do not match" )
  
  call delete_file( file_name )
end subroutine

!> Unit tests for [[occupy_cdft]]
subroutine test_occupy_cdft( test_report )
  !> Unit test report
  type(unit_test_type) :: test_report

  real(dp), parameter :: tol = 1e-8_dp
  integer(i32), parameter :: n_exc = 3
  integer(i32), parameter :: n_kpt = 2
  integer(i32), parameter :: n_valence_states = 2
  integer(i32), parameter :: n_conduction_states = 4
  integer(i32), parameter :: n_states = n_valence_states + n_conduction_states
  logical :: spin_polarization
  integer(i32) :: i, iv(n_exc), ic(n_exc), ik(n_exc)
  complex(dp), parameter :: z(n_exc) = [ (0.17_dp, 0.26_dp), (0.4_dp, 0.0_dp), (0.0_dp, 0.7_dp) ]
  real(dp) :: occupation_factors(n_states, n_kpt), occupation_factors_ref(n_states, n_kpt), wkpt(n_kpt)
  real(dp) :: occupation_factors_spin(2*n_states, n_kpt), occupation_factors_spin_ref(2*n_states, n_kpt)
  type(ExcitonCoefficients) :: exc

  ! 1st test (spin unpolarized)
  spin_polarization = .false.  
  wkpt = spread( 1.0_dp/n_kpt, dim=1, ncopies=n_kpt )
  occupation_factors_ref = 0._dp
  occupation_factors_ref(1:n_valence_states, 1:n_kpt) = 2._dp
  occupation_factors = occupation_factors_ref
  ik = [1, 2, 1]
  iv = [1, 2, 2]
  ic = [3, 4, 6]
  call exc%set_attributes( ik, iv, ic, z )
  do i = 1, n_exc
    call change_occupation( iv(i), ic(i), abs(z(i))**2/wkpt(ik(i)), occupation_factors_ref(:, ik(i)) )
  end do
  call occupy_cdft( exc, wkpt, occupation_factors )
  call test_report%assert( all_close(occupation_factors, occupation_factors_ref, tol), "Wrong occupation factors in test 1" )

  ! 2nd test (spin polarized)
  spin_polarization = .true.
  occupation_factors_spin_ref = 0._dp
  occupation_factors_spin_ref(1:n_valence_states, 1:n_kpt) = 1._dp
  occupation_factors_spin_ref((n_states+1):(n_states+n_valence_states), 1:n_kpt) = 1._dp
  occupation_factors_spin = occupation_factors_spin_ref
  ik = [1, 2, 1]
  iv = [1, 2, 8]
  ic = [3, 4, 12]
  call exc%set_attributes( ik, iv, ic, z )
  do i = 1, n_exc
    call change_occupation( iv(i), ic(i), abs(z(i))**2/wkpt(ik(i)), occupation_factors_spin_ref(:, ik(i)) )
  end do
  call occupy_cdft( exc, wkpt, occupation_factors_spin )
  call test_report%assert( all_close(occupation_factors_spin, occupation_factors_spin_ref, tol), "Wrong occupation factors in test 2" )

  contains 
    pure subroutine change_occupation( i_sub, i_add, delta, occupation )
      integer(i32), intent(in) :: i_sub, i_add
      real(dp), intent(in) :: delta
      real(dp), intent(inout) :: occupation(:)

      occupation(i_sub) = occupation(i_sub) - delta
      occupation(i_add) = occupation(i_add) + delta
    end subroutine

end subroutine

!> Delete a file, if it exists
subroutine delete_file( file_name )
  !> File name to delete
  character(len=*), intent(in) :: file_name

  integer(i32) :: unit
  logical :: file_exists

  ! Delete test file
  inquire( file=trim(file_name), exist=file_exists )
  if( file_exists ) then
    open( newunit=unit, file=trim(file_name), status='old' )
    close( unit, status='delete' )
  end if
  
end subroutine

end module