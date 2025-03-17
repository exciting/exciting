module rttddft_NumberExcitations_test
  use constants, only: zone, zi, sqrt_two
  use exciting_mpi, only: mpiinfo
  use math_utils, only: all_close
  use mock_arrays, only: complex_positive_definite_matrix_5x5, complex_hermitian_matrix_5x5
  use precision, only: dp, i32
  use rttddft_NumberExcitations, only: obtain_number_excitations
  use unit_test_framework, only : unit_test_type
  use xlapack, only: solve_generalized_hermitian_eigenproblem
  use rttddft_Wavefunction, only: wavefunction_set, initialize_wavefunction_set

  implicit none

  private

  public :: rttddft_NumberExcitations_test_driver

  real(dp), parameter :: tol = 1.0e-10_dp
 
contains

  subroutine rttddft_NumberExcitations_test_driver( mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional, intent(in) :: kill_on_failure
    
    type(unit_test_type) :: test_report
    integer(i32), parameter :: n_assertions_test_obtain_number_excitations = 6
    integer(i32), parameter :: n_assertions = n_assertions_test_obtain_number_excitations

    character(len=*), parameter :: module_tested = 'rttddft_NumberExcitations'

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_obtain_number_excitations( mpiglobal, test_report )

    ! report results
    if ( present( kill_on_failure ) ) then
      call test_report%report( module_tested, kill_on_failure )
    else
      call test_report%report( module_tested )
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine

  !> Unit tests for [[obtain_number_excitations]]
  subroutine test_obtain_number_excitations( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_identifier = "test_obtain_number_excitations"
    integer(i32), parameter :: n_states_gnd = 5
    integer(i32), parameter :: n_states = 3
    integer(i32), parameter :: n_frozen = 1
    integer(i32), parameter :: i_VBM = 2
    integer(i32), parameter :: i_CBm = i_VBM + 1
    integer(i32), parameter :: n_dim = size( complex_positive_definite_matrix_5x5, 1 )
    integer(i32) :: ik, n_kpt, n_kpt_per_rank, first_k, last_k
    real(dp) :: n_exc, n_gs, n_exc_ref, n_gs_ref, eigs_gnd(n_states_gnd), eigs(n_states)
    real(dp), allocatable :: wkpt(:), occ_gnd(:, :), occ(:)
    complex(dp), allocatable :: H_gnd(:, :), H(:, :), S(:, :), overlap(:, :, :)
    complex(dp), allocatable :: psi_gnd(:, :, :), tmp(:, :), proj(:, :), psi_t(:, :, :)
    class(wavefunction_set), allocatable :: psi, psi_no_frozen

    ! Initialization
    n_kpt_per_rank = 2
    n_kpt = n_kpt_per_rank*( mpiglobal%procs )
    wkpt = [ (ik**2, ik = 1, n_kpt) ]
    wkpt = wkpt/sum( wkpt )
    first_k = n_kpt_per_rank*( mpiglobal%rank ) + 1
    last_k = first_k + n_kpt_per_rank - 1
    allocate( overlap(n_dim, n_dim, n_kpt), occ_gnd(n_states_gnd, n_kpt) )
    allocate( psi_gnd(n_dim, n_states_gnd, n_kpt), psi_t(n_dim, n_frozen + n_states, n_kpt) )
    H_gnd = complex_hermitian_matrix_5x5
    do ik = 1, n_kpt
      overlap(:, :, ik) = complex_positive_definite_matrix_5x5
      S = overlap(:, :, ik); H = H_gnd
      call solve_generalized_hermitian_eigenproblem( H, S, tol, eigs_gnd, psi_gnd(:, :, ik) )
      S = overlap(:, :, ik); H = H_gnd + ik*S
      call solve_generalized_hermitian_eigenproblem( H, S, tol, eigs, psi_t(:, :, ik) )
      occ_gnd(1:i_VBM, ik) = 2._dp 
      occ_gnd(i_CBm:, ik) = 0._dp
      if( modulo(n_kpt, ik) == 1 ) then
        occ_gnd(i_VBM, ik) = occ_gnd(i_VBM, ik) - real(ik, dp)/n_kpt
        occ_gnd(i_CBm, ik) = occ_gnd(i_CBm, ik) + real(ik, dp)/n_kpt
      end if
    end do
    call initialize_wavefunction_set( psi, .true., .false., n_frozen, psi_gnd, occ_gnd, tol )
    psi%active = psi_t(:, n_frozen + 1:, :)

    ! Obtain ref
    n_exc_ref = 0._dp; n_gs_ref = 0._dp
    do ik = 1, n_kpt
      tmp = matmul( overlap(:, :, ik), psi%groundstate(:, :, ik) )
      proj = matmul( conjg( transpose( tmp ) ), psi%active(:, :, ik) )
      occ = matmul( abs((proj(:, :))**2), occ_gnd(n_frozen + 1 : n_frozen + n_states, ik) )
      occ(1 : n_frozen) = occ(1 : n_frozen) + occ_gnd(1 : n_frozen, ik)
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      n_exc_ref = n_exc_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) <= tol )
    end do
    ! Obtain n_exc and n_gs
    call obtain_number_excitations( psi, overlap(:, :, first_k:last_k), tol, &
      occ_gnd(:, first_k:last_k), wkpt(first_k:last_k), mpiglobal, n_exc, n_gs )
    
    call test_report%assert( all_close( n_exc, n_exc_ref, tol ), test_identifier//' - calculated n_exc does not match reference.' )
    call test_report%assert( all_close( n_gs, n_gs_ref, tol ), test_identifier//' - calculated n_gs does not match reference.' )

    ! Test case with no excited states
    psi%active = psi%groundstate(:, psi%first_active(): psi%n_occupied(), :)
    call obtain_number_excitations( psi, overlap(:, :, first_k:last_k), tol, &
      occ_gnd(:, first_k:last_k), wkpt(first_k:last_k), mpiglobal, n_exc, n_gs )
    
    call test_report%assert( all_close( n_exc, 0._dp, tol ), test_identifier//' - calculated n_exc does not match reference.' )
    call test_report%assert( all_close( n_gs, sum(occ_gnd(:, 1)), tol ), test_identifier//' - calculated n_gs does not match reference.' )

    ! Test case with no frozen states
    call initialize_wavefunction_set( psi_no_frozen, .true., .false., 0, psi_gnd, occ_gnd, tol )
    psi_no_frozen%active = psi_t(:, 1 : n_states, :)

    ! Obtain ref
    n_exc_ref = 0._dp; n_gs_ref = 0._dp
    do ik = 1, n_kpt
      tmp = matmul( overlap(:, :, ik), psi_no_frozen%groundstate(:, :, ik) )
      proj = matmul( conjg( transpose( tmp ) ), psi_no_frozen%active(:, :, ik) )
      occ = matmul( abs((proj(:, :))**2), occ_gnd(1 : n_states, ik) )
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      n_exc_ref = n_exc_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) <= tol )
    end do
    ! Obtain n_exc and n_gs
    call obtain_number_excitations( psi_no_frozen, overlap(:, :, first_k:last_k), tol, &
      occ_gnd(:, first_k:last_k), wkpt(first_k:last_k), mpiglobal, n_exc, n_gs )
    
    call test_report%assert( all_close( n_exc, n_exc_ref, tol ), test_identifier//' - calculated n_exc does not match reference.' )
    call test_report%assert( all_close( n_gs, n_gs_ref, tol ), test_identifier//' - calculated n_gs does not match reference.' )

  end subroutine


end module