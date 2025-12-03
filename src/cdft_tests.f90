!> Module that contains the CDFT unit tests
module cdft_tests
  use constants, only: zi, zone, zzero
  use cdft, only: cdft_input_keys, &
                  deallocate_cdft_global_arrays, &
                  determine_cdft_occupations, &
                  ExcitonCoefficients, &
                  initialize_cdft_global_arrays, &
                  Occupations, &
                  set_overlap_times_psi_gs, &
                  update_occupations_with_the_maximum_overlap_method
  use file_utils, only: delete_file
  use math_utils, only: all_close
  use modmpi, only: mpiinfo
  use mock_arrays, only: complex_unitary_matrix_5x5
  use precision, only: i32, dp
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none

  private

  public :: run_cdft_test_driver

  character(len=*), parameter :: warnings = "WARNINGS.OUT"

contains 
!> Run the CDFT unit tests
subroutine run_cdft_test_driver( mpiglobal, kill_on_failure )
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal
  !> Kill the program before the test driver finishes if an assertion fails
  logical, optional :: kill_on_failure
  
  type(unit_test_type) :: test_report
  character(len=*), parameter :: test_driver_name = "cdft"

  call test_report%init( mpiglobal )

  ! Run and assert tests
  call test_ExcitonCoefficients_get_from_file( test_report, mpiglobal )
  call test_Occupations_get_from_file( test_report, mpiglobal )
  call test_determine_cdft_occupations( test_report, mpiglobal )
  call test_occupy_update_occupations_max_overl_meth( test_report )
  
  call test_report%report( test_driver_name, kill_on_failure )
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
  integer(i32) :: i, unit, i_err
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
  
  call delete_file( file_name, i_err )
  if( mpiglobal%is_root ) call delete_file( warnings, i_err )
end subroutine

!> Unit tests for [[cdft(module):Occupations_get_from_file(subroutine)]]
subroutine test_Occupations_get_from_file( test_report, mpiglobal )
  !> Unit test report
  type(unit_test_type) :: test_report
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal

  character(len=*), parameter :: test_name = "test_Occupations_get_from_file: "
  integer(i32), parameter :: n_non_zero_ref = 4
  integer(i32), parameter :: idx_states_ref(n_non_zero_ref) = [4, 3, 8, 9] 
  integer(i32), parameter :: idx_kpt_ref(n_non_zero_ref) = [1, 1, 7, 9] 
  character(len=*), parameter :: fake_name = "test"
  real(dp), parameter :: occ_ref(n_non_zero_ref) = [ 1.8_dp, 1.7_dp, 0.4_dp, 0.1_dp]
  real(dp), parameter :: tol = 1e-8_dp
  integer(i32) :: i, unit, i_err
  character(len=:), allocatable :: file_name
  integer(i32), allocatable :: idx_states(:), idx_kpt(:)
  real(dp), allocatable :: occ_factors(:)
  type(Occupations) :: occ

  ! Each MPI rank reads/writes its own file
  file_name = fake_name // to_char( mpiglobal%rank )
  ! Write to file
  open( newunit = unit, file = trim(file_name), action = 'write' )
  write( unit, * ) n_non_zero_ref
  do i = 1, n_non_zero_ref
    write( unit, '(2I6, F18.10)' ) idx_states_ref(i), idx_kpt_ref(i), occ_ref(i)
  end do
  close( unit )

  call occ%get_from_file( file_name )
  call occ%get_attributes( idx_kpt, idx_states, occ_factors )
  call test_report%assert( n_non_zero_ref == size( occ_factors ), test_name // "Wrong number of occupation factors" )
  call test_report%assert( all( idx_kpt == idx_kpt_ref ), test_name // "k-point indexes do not match" )
  call test_report%assert( all( idx_states == idx_states_ref ), test_name // "state indexes do not match" )
  call test_report%assert( all_close(occ_factors, occ_ref, tol), test_name // "occupation factors do not match" )

  call delete_file( file_name, i_err )
  if( mpiglobal%is_root ) call delete_file( warnings, i_err )
end subroutine

!> Unit tests for [[cdft(module):determine_cdft_occupations(subroutine)]]
subroutine test_determine_cdft_occupations( test_report, mpiglobal )
  !> Unit test report
  type(unit_test_type) :: test_report
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal

  real(dp), parameter :: tol = 1e-8_dp
  integer(i32), parameter :: n_exc = 3
  integer(i32), parameter :: n_kpt = 2
  integer(i32), parameter :: n_valence_states = 2
  integer(i32), parameter :: n_conduction_states = 4
  integer(i32), parameter :: n_states = n_valence_states + n_conduction_states
  logical :: spin_polarization
  integer(i32) :: i, iv(n_exc), ic(n_exc), ik(n_exc), is(n_exc), i_err
  complex(dp), parameter :: z(n_exc) = [ (0.17_dp, 0.26_dp), (0.4_dp, 0.0_dp), (0.0_dp, 0.7_dp) ]
  real(dp) :: occupation_changed(n_exc)
  real(dp) :: occupation_factors(n_states, n_kpt), occupation_factors_ref(n_states, n_kpt), wkpt(n_kpt)
  real(dp) :: occupation_factors_spin(2*n_states, n_kpt), occupation_factors_spin_ref(2*n_states, n_kpt)
  type(cdft_input_keys) :: cdft_inp
  type(ExcitonCoefficients) :: exc
  type(Occupations) :: occ

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
  call cdft_inp%mock( exc )
  call determine_cdft_occupations( cdft_inp, wkpt, occupation_factors )
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
  call cdft_inp%mock( exc )
  call determine_cdft_occupations( cdft_inp, wkpt, occupation_factors_spin )
  call test_report%assert( all_close(occupation_factors_spin, occupation_factors_spin_ref, tol), "Wrong occupation factors in test 2" )

  ! 3rd test (Occupations)
  occupation_factors_ref = 0._dp
  occupation_factors_ref(1:n_valence_states, 1:n_kpt) = 2._dp
  occupation_factors = occupation_factors_ref
  ik = [1, 2, 1]
  is = [1, 2, 5]
  occupation_changed = [1.5_dp, 1.7_dp, 0.8_dp]
  call occ%set_attributes( ik, is, occupation_changed )
  call cdft_inp%mock( occ )
  call determine_cdft_occupations( cdft_inp, wkpt, occupation_factors )
  do i = 1, n_exc
    occupation_factors_ref(is(i), ik(i)) = occupation_changed(i)
  end do
  call test_report%assert( all_close(occupation_factors, occupation_factors_ref, tol), "Wrong occupation factors in test 3" )
  if( mpiglobal%is_root ) call delete_file( warnings, i_err )

  contains 
    pure subroutine change_occupation( i_sub, i_add, delta, occupation )
      integer(i32), intent(in) :: i_sub, i_add
      real(dp), intent(in) :: delta
      real(dp), intent(inout) :: occupation(:)

      occupation(i_sub) = occupation(i_sub) - delta
      occupation(i_add) = occupation(i_add) + delta
    end subroutine

end subroutine

!> Unit tests for [[update_occupations_with_the_maximum_overlap_method]]
subroutine test_occupy_update_occupations_max_overl_meth( test_report )
  !> Unit test report
  type(unit_test_type) :: test_report
  real(dp), parameter :: tol = 1e-8_dp
  integer(i32), parameter :: n_kpt = 1
  integer(i32), parameter :: first_kpt = 1
  integer(i32), parameter :: n_basis = size( complex_unitary_matrix_5x5, 1 )
  integer(i32), parameter :: n_val = 3
  integer(i32), parameter :: n_states = size( complex_unitary_matrix_5x5, 2 )
  integer(i32) :: i, ik
  integer(i32), allocatable :: indexes(:)
  real(dp), allocatable :: occ_gs(: , :), occ_expected(:, :)
  complex(dp), allocatable :: psi_gs(:, :, :), psi(:, :, :), S(:, :)
  
  ! Initialization
  allocate( psi_gs(n_basis, n_states, n_kpt), occ_gs(n_states, n_kpt) )
  do ik = 1, n_kpt
    psi_gs(:, :, ik) = complex_unitary_matrix_5x5
    occ_gs(:, ik) = [ ( merge(2._dp, 0._dp, i<=n_val), i = 1, n_states) ]
  end do

  ! Identity matrix
  allocate( S(n_basis, n_basis), source=zzero )
  do i = 1, n_basis
    S(i, i) = zone
  end do

  ! Exchange states
  allocate( indexes(n_states), psi(n_basis, n_states, n_kpt), occ_expected(n_states, n_kpt) )
  do ik = 1, n_kpt
    indexes = [ (modulo(n_states - i + ik, n_states) + 1, i = 1, n_states) ]
    occ_expected(:, ik) = occ_gs(indexes, ik)
    psi(:, :, ik) = psi_gs(:, indexes, ik)
  end do

  call initialize_cdft_global_arrays( psi_gs, first_kpt )
  do ik = 1, n_kpt
    call set_overlap_times_psi_gs( ik, S )
  end do
  call update_occupations_with_the_maximum_overlap_method( psi, occ_gs )
  call deallocate_cdft_global_arrays( )
  call test_report%assert( all_close(occ_expected, occ_gs, tol), &
    "cdft: occupations do not match; max. diff = " // to_char( maxval( abs(occ_expected-occ_gs) ) ) )
end subroutine

end module