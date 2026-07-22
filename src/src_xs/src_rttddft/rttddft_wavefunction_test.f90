module rttddft_Wavefunction_test
  use constants, only: zzero,zone, zi, sqrt_two
  use exciting_mpi, only: mpiinfo
  use math_utils, only: all_close
  use mock_arrays, only: complex_hermitian_matrix_5x5, complex_matrix_7x5, &
    complex_positive_definite_matrix_5x5, complex_unitary_matrix_5x5
  use mod_kpointset, only: k_set
  use normalize, only: norm_squared_with_positive_matrix
  use precision, only: dp, i32
  use rttddft_Overlap, only: overlap_set
  use rttddft_Wavefunction, only: obtain_occupations, wavefunction_set, initialize_wavefunction_set
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type
  use xlapack, only: solve_generalized_hermitian_eigenproblem, norm

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

    character(len=*), parameter :: module_tested = 'rttddft_Wavefunction'

    ! Initialize test object
    call test_report%init( mpiglobal )

    ! Run and assert tests
    call test_obtain_occupations( test_report )
    call test_obtain_number_excitations( mpiglobal, test_report )

    ! report results
    call test_report%report( module_tested, kill_on_failure )

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
    integer(i32), parameter :: n_states = 3 ! occupied states
    integer(i32), parameter :: n_frozen = 1
    integer(i32), parameter :: i_VBM = 2 
    integer(i32), parameter :: i_CBm = i_VBM + 1 ! metallic system: iCBm = n_states
    integer(i32), parameter :: n_dim = 7
    integer(i32) :: ik, j, n_kpt, n_kpt_per_rank, first_k, last_k, test_counter
    real(dp) :: n_exc, n_gs, n_t_ref, n_exc_ref, n_gs_ref
    real(dp) :: eigs_gnd(n_states_gnd), eigs(n_frozen + n_states)
    real(dp), allocatable :: wkpt(:), occ_gnd(:, :), occ(:), norms_squared(:)
    complex(dp), allocatable :: H_gnd(:, :), H(:, :), S_aux(:, :), aux(:, :), overlap(:, :, :)
    complex(dp), allocatable :: psi_gnd(:, :, :), tmp(:, :), proj(:, :), psi_t(:, :, :)
    type(overlap_set) :: S
    class(wavefunction_set), allocatable :: psi, psi_no_frozen
    type(k_set) :: kset

    ! Initialization: k-points
    n_kpt_per_rank = 2
    n_kpt = n_kpt_per_rank*( mpiglobal%procs )
    wkpt = [ (ik**2, ik = 1, n_kpt) ]
    wkpt = wkpt/sum( wkpt )
    kset%wkpt = wkpt
    first_k = n_kpt_per_rank*( mpiglobal%rank ) + 1
    last_k = first_k + n_kpt_per_rank - 1
    ! Initialization: H and S
    allocate( occ_gnd(n_states_gnd, n_kpt), overlap(n_dim, n_dim, n_kpt) )
    allocate( psi_gnd(n_dim, n_states_gnd, n_kpt), psi_t(n_dim, n_frozen + n_states, n_kpt) )
    allocate( H_gnd(n_dim, n_dim), H(n_dim, n_dim), S_aux(n_dim, n_dim), aux(n_dim, n_dim) )
    H_gnd = zzero
    H_gnd(1:5, 1:5) = complex_hermitian_matrix_5x5
    aux = zzero
    aux(1:7, 1:5) = complex_matrix_7x5
    aux(1:7, 6:7) = complex_matrix_7x5(1:7, 1:2)
    S_aux = transpose( aux ) ! auxiliary operation
    aux = conjg( S_aux ) + aux
    S_aux = zzero
    S_aux(1:5, 1:5) = complex_positive_definite_matrix_5x5
    S_aux(6, 6) = 1._dp; S_aux(7, 7) = 856.71000329_dp
    do ik = 1, n_kpt
      overlap(:, :, ik) = S_aux
    end do
    call S%allocate( .true., n_dim, first_k, last_k )
    S%array = overlap(:, :, first_k:last_k)
    ! Initialization: wavefunctions and occupations
    do ik = 1, n_kpt
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
    ! Change the norm of one wavefunction
    psi_t(:, i_CBm, 1) = 0.97321_dp*psi_t(:, i_CBm, 1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test case: LAPWlo basis with frozen states
    call initialize_wavefunction_set( psi, .true., first_k, kset, .false., n_frozen, &
      psi_gnd(:, :, first_k:last_k), occ_gnd(:, first_k:last_k), tol )
    psi%active(:, :, first_k:last_k) = psi_t(:, n_frozen + 1:i_CBM, first_k:last_k)

    ! Obtain ref
    n_exc_ref = 0._dp; n_gs_ref = 0._dp; n_t_ref = 0._dp
    allocate( norms_squared(n_frozen + n_states), source = 1._dp )
    do ik = 1, n_kpt
      tmp = conjg( matmul( overlap(:, :, ik), psi_gnd(:, :, ik) ) )
      proj = matmul( transpose( tmp ), psi_t(:, :, ik) )
      occ = matmul( abs((proj(:, n_frozen + 1 : n_frozen + n_states))**2), &
        occ_gnd(n_frozen + 1 : n_frozen + n_states, ik) )
      occ(1 : n_frozen) = occ(1 : n_frozen) + occ_gnd(1 : n_frozen, ik)
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      call norm_squared_with_positive_matrix( psi_t(:, :, ik), overlap(:, :, ik), norms_squared )
      n_t_ref = n_t_ref + wkpt(ik)*dot_product( occ_gnd(1 : n_frozen + n_states, ik), norms_squared )
    end do
    n_exc_ref = n_t_ref - n_gs_ref
    ! Obtain n_exc and n_gs
    test_counter = 1
    call psi%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )
    call test_asserts( test_counter )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test case: LAPWlo basis with frozen states and no excited states
    test_counter = test_counter + 1
    psi%active(:, :, :) = psi%groundstate_lapwlo(:, psi%first_active(): psi%n_occupied(), :)
    call psi%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )
    n_exc_ref = 0._dp; n_gs_ref = sum(occ_gnd(:, 1))
    call test_asserts( test_counter )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test case: LAPWlo basis with no frozen states
    test_counter = test_counter + 1
    call initialize_wavefunction_set( psi_no_frozen, .true., first_k, kset, .false., &
      0, psi_gnd(:, :, first_k:last_k), occ_gnd(:, first_k:last_k), tol )
    psi_no_frozen%active(:, 1 : n_states, first_k:last_k) = psi_t(:, 1 : n_states, first_k:last_k)

    ! Obtain ref
    n_t_ref = 0._dp; n_gs_ref = 0._dp
    deallocate( proj, occ, tmp )
    do ik = 1, n_kpt
      tmp = matmul( overlap(:, :, ik), psi_gnd(:, :, ik) )
      proj = matmul( conjg( transpose( tmp ) ), psi_t(:, 1 : n_states, ik) )
      occ = matmul( abs((proj(:, :))**2), occ_gnd(1 : n_states, ik) )
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      call norm_squared_with_positive_matrix( psi_t(:, :, ik), overlap(:, :, ik), norms_squared )
      n_t_ref = n_t_ref + wkpt(ik)*dot_product( occ_gnd(1:n_states, ik), norms_squared(1:n_states) )
    end do
    n_exc_ref = n_t_ref - n_gs_ref
    ! Obtain n_exc and n_gs
    call psi_no_frozen%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )
    call test_asserts( test_counter )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test case: KS basis with frozen states
    test_counter = test_counter + 1
    deallocate( psi )
    call initialize_wavefunction_set( psi, .false., first_k, kset, .false., n_frozen, &
      psi_gnd(:, :, first_k:last_k), occ_gnd(:, first_k:last_k), tol )
    do ik = first_k, last_k
      psi%active(:, :, ik) = complex_unitary_matrix_5x5(1:5, 1:n_states-n_frozen) * &
        cmplx( cos( 1.0_dp*ik ), sin( 1.0_dp*ik ), dp )
    end do
    ! Obtain ref
    n_t_ref = 0._dp; n_gs_ref = 0._dp
    deallocate( proj, occ )
    do ik = 1, n_kpt
      proj = complex_unitary_matrix_5x5(1:5, 1:n_states-n_frozen)*cmplx( cos( 1.0_dp*ik ), sin( 1.0_dp*ik ), dp )
      occ = matmul( abs((proj(:, :))**2), occ_gnd(n_frozen + 1 : n_states, ik) )
      occ(1 : n_frozen) = occ(1 : n_frozen) + occ_gnd(1 : n_frozen, ik)
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      norms_squared = 1._dp
      do j = 1, n_states - n_frozen
        norms_squared(j+n_frozen) = norm( proj(:, j) )**2
      end do
      n_t_ref = n_t_ref + wkpt(ik)*dot_product( occ_gnd(1:n_frozen+n_states, ik), norms_squared )
    end do
    n_exc_ref = n_t_ref - n_gs_ref
    ! Obtain n_exc and n_gs
    call psi%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )
    call test_asserts( test_counter )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test case: KS basis with no frozen states
    test_counter = test_counter + 1
    deallocate( psi_no_frozen )
    call initialize_wavefunction_set( psi_no_frozen, .false., first_k, kset, .false., &
      0, psi_gnd(:, :, first_k:last_k), occ_gnd(:, first_k:last_k), tol )
    do ik = first_k, last_k
      psi_no_frozen%active(:, :, ik) = complex_unitary_matrix_5x5(1:5, 1:n_states) * &
        cmplx( cos( 1.0_dp*ik ), sin( 1.0_dp*ik ), dp )
    end do
    ! Obtain ref
    n_t_ref = 0._dp; n_gs_ref = 0._dp
    deallocate( proj, occ )
    do ik = 1, n_kpt
      proj = complex_unitary_matrix_5x5(1:5, 1:n_states)*cmplx( cos( 1.0_dp*ik ), sin( 1.0_dp*ik ), dp )
      occ = matmul( abs((proj(:, :))**2), occ_gnd(1:n_states, ik) )
      n_gs_ref = n_gs_ref + wkpt(ik)*sum( occ, occ_gnd(:, ik) >= tol )
      do j = 1, n_states
        norms_squared(j) = norm( proj(:, j) )**2
      end do
      n_t_ref = n_t_ref + wkpt(ik)*dot_product( occ_gnd(1:n_states, ik), norms_squared(1:n_states) )
    end do
    n_exc_ref = n_t_ref - n_gs_ref
    ! Obtain n_exc and n_gs
    call psi_no_frozen%obtain_number_excitations( S, mpiglobal, n_exc, n_gs )
    call test_asserts( test_counter )

    contains
      !> Auxiliary subroutine to call assertions
      subroutine test_asserts( counter )
        integer(i32), intent(in) :: counter
        call test_assert_all_close( 'n_exc', n_exc, n_exc_ref, counter )
        call test_assert_all_close( 'n_gs', n_gs, n_gs_ref, counter )
      end subroutine

      !> Auxiliary subroutine to assert that calculated and reference values are close
      subroutine test_assert_all_close( name, n, n_ref, counter )
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: n, n_ref
        integer(i32), intent(in) :: counter
        call test_report%assert( all_close( n, n_ref, tol ), report_message( test_identifier, name, counter ) )
      end subroutine
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
