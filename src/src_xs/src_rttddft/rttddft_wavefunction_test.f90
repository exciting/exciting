module rttddft_Wavefunction_test
  use constants, only: zone, zi, sqrt_two
  use exciting_mpi, only: mpiinfo
  use math_utils, only: all_close
  use mock_arrays, only: complex_positive_definite_matrix_5x5, complex_hermitian_matrix_5x5, complex_matrix_5x5
  use mod_kpointset, only: k_set
  use precision, only: dp, i32
  use rttddft_Overlap, only: overlap_set
  use rttddft_Wavefunction, only: obtain_occupations, wavefunction_set, initialize_wavefunction_set
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type
  use xlapack, only: solve_generalized_hermitian_eigenproblem

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
    integer(i32), parameter :: n_assertions_test_obtain_occupations = 8
    integer(i32), parameter :: n_assertions = n_assertions_test_obtain_occupations

    character(len=*), parameter :: module_tested = 'rttddft_Wavefunction'

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_obtain_occupations( test_report )
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

  
  subroutine test_obtain_occupations( test_report )
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

  !> Unit tests for [[obtain_number_excitations]]
  subroutine test_obtain_number_excitations( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_identifier = "test_obtain_number_excitations"
    integer(i32), parameter :: n_states_gnd = 5
    integer(i32), parameter :: n_states = 3 ! occupied + empty states
    integer(i32), parameter :: n_frozen = 1
    integer(i32), parameter :: i_VBM = 2 
    integer(i32), parameter :: i_CBm = i_VBM + 1 
    integer(i32), parameter :: n_dim = size( complex_positive_definite_matrix_5x5, 1 )
    integer(i32) :: ik, n_kpt, n_kpt_per_rank, first_k, last_k, test_counter
    real(dp) :: n_exc, n_gs, n_exc_ref, n_gs_ref, eigs_gnd(n_states_gnd), eigs(n_frozen + n_states)
    real(dp), allocatable :: wkpt(:), occ_gnd(:, :), occ(:)
    complex(dp), allocatable :: H_gnd(:, :), H(:, :), S_aux(:, :), aux(:, :), overlap(:, :, :)
    type(overlap_set) :: S
    complex(dp), allocatable :: psi_gnd(:, :, :), tmp(:, :), proj(:, :), psi_t(:, :, :)
    class(wavefunction_set), allocatable :: psi, psi_no_frozen
    type(k_set) :: kset

    ! Initialization
    n_kpt_per_rank = 2
    n_kpt = n_kpt_per_rank*( mpiglobal%procs )
    wkpt = [ (ik**2, ik = 1, n_kpt) ]
    wkpt = wkpt/sum( wkpt )
    kset%wkpt = wkpt
    first_k = n_kpt_per_rank*( mpiglobal%rank ) + 1
    last_k = first_k + n_kpt_per_rank - 1
    allocate( occ_gnd(n_states_gnd, n_kpt), overlap(n_dim, n_dim, n_kpt) )
    allocate( psi_gnd(n_dim, n_states_gnd, n_kpt), psi_t(n_dim, n_frozen + n_states, n_kpt) )
    H_gnd = complex_hermitian_matrix_5x5
    aux = conjg(complex_matrix_5x5)
    aux = transpose(aux) + complex_matrix_5x5
    do ik = 1, n_kpt
      overlap(:, :, ik) = complex_positive_definite_matrix_5x5
      S_aux = overlap(:, :, ik); H = H_gnd
      call solve_generalized_hermitian_eigenproblem( H, S_aux, tol, eigs_gnd, psi_gnd(:, :, ik) )
      S_aux = overlap(:, :, ik); H = H_gnd + ik*aux
      call solve_generalized_hermitian_eigenproblem( H, S_aux, tol, eigs, psi_t(:, :, ik) )
      occ_gnd(1:i_VBM, ik) = 2._dp
      occ_gnd(i_CBm:, ik) = 0._dp
    end do
    !!! Metallic system: i_CBM will be occupied
    occ_gnd(i_VBM, 1) = occ_gnd(i_VBM, 1) - 1.0_dp/n_kpt
    occ_gnd(i_CBm, 1) = occ_gnd(i_CBm, 1) + 1.0_dp/n_kpt
    call S%allocate( .true., n_dim, first_k, last_k )
    S%array = overlap(:, :, first_k:last_k)
    call initialize_wavefunction_set( psi, .true., first_k, kset, .false., n_frozen, psi_gnd(:, :, first_k:last_k), occ_gnd(:, first_k:last_k), tol )
    psi%active(:, :, first_k:last_k) = psi_t(:, n_frozen + 1:i_CBM, first_k:last_k)

    ! Obtain ref
    n_exc_ref = 0._dp; n_gs_ref = 0._dp
    do ik = 1, n_kpt
      tmp = matmul( overlap(:, :, ik), psi_gnd(:, :, ik) )
      proj = matmul( conjg( transpose( tmp ) ), psi_t(:, :, ik) )
      occ = matmul( abs((proj(:, n_frozen + 1 : n_frozen + n_states))**2), occ_gnd(n_frozen + 1 : n_frozen + n_states, ik) )
      occ(1 : n_frozen) = occ(1 : n_frozen) + occ_gnd(1 : n_frozen, ik)
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      n_exc_ref = n_exc_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) <= tol )
    end do
    ! Obtain n_exc and n_gs
    test_counter = 1
    call psi%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )
    call test_report%assert( all_close( n_exc, n_exc_ref, tol ), report_message( test_identifier, 'n_exc', test_counter ) )
    call test_report%assert( all_close( n_gs, n_gs_ref, tol ), report_message( test_identifier, 'n_gs', test_counter ) )

    ! Test case: frozen states and no excited states
    test_counter = test_counter + 1
    psi%active(:, :, :) = psi%groundstate(:, psi%first_active(): psi%n_occupied(), :)
    call psi%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )

    call test_report%assert( all_close( n_exc, 0._dp, tol ), report_message( test_identifier, 'n_exc', test_counter ) )
    call test_report%assert( all_close( n_gs, sum(occ_gnd(:, 1)), tol ), report_message( test_identifier, 'n_gs', test_counter ) )

    ! Test case with no frozen states
    test_counter = test_counter + 1
    call initialize_wavefunction_set( psi_no_frozen, .true., first_k, kset, .false., 0, psi_gnd(:, :, first_k:last_k), occ_gnd(:, first_k:last_k), tol )
    psi_no_frozen%active(:, 1 : n_states, first_k:last_k) = psi_t(:, 1 : n_states, first_k:last_k)

    ! Obtain ref
    n_exc_ref = 0._dp; n_gs_ref = 0._dp
    deallocate( proj, occ, tmp )
    do ik = 1, n_kpt
      tmp = matmul( overlap(:, :, ik), psi_gnd(:, :, ik) )
      proj = matmul( conjg( transpose( tmp ) ), psi_t(:, 1 : n_states, ik) )
      occ = matmul( abs((proj(:, :))**2), occ_gnd(1 : n_states, ik) )
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      n_exc_ref = n_exc_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) <= tol )
    end do
    ! Obtain n_exc and n_gs
    call psi_no_frozen%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )

    call test_report%assert( all_close( n_exc, n_exc_ref, tol ), report_message( test_identifier, 'n_exc', test_counter ) )
    call test_report%assert( all_close( n_gs, n_gs_ref, tol ), report_message( test_identifier, 'n_gs', test_counter ) )
  end subroutine

  !> (private) generates a report message, given a test case and test number
  function report_message( test_id, test_case, test_number ) result(message)
    !> Name of the subroutine calling `report_message`
    character(len=*), intent(in) :: test_id
    !> Test (sub)identifier, pointing out which case in [[test_id]] is tested
    character(len=*), intent(in) :: test_case
    !> Number for the test
    integer(i32), intent(in) :: test_number
    !> Message to return
    character(len=:), allocatable :: message
    message = test_id // ' - ' // test_case // ' does not match reference: test - ' // to_char(test_number)
  end function

end module
