module propagators_test
  use constants, only: zone, zzero, zi
  use exciting_mpi, only: mpiinfo
  use integration, only: rk4 => ODESolver_RungeKutta4thOrder
  use math_utils, only: all_close
  use matrix_exp, only: &
    & exp_hermitian => exp_hermitian_matrix_times_vectors, & 
    & exp_general => exp_general_matrix_times_vectors, &
    & exp_houston => exphouston_hermitian_matrix_times_vectors
  use mock_arrays, only: herm => complex_hermitian_matrix_5x5, &
                         pos => complex_positive_definite_matrix_5x5, &
                         c5x7 => complex_matrix_5x7, &
                         c5x5 => complex_matrix_5x5, &
                         c7x5 => complex_matrix_7x5
  use precision, only: dp, i32
  use propagators, only: propagator, create_propagator, &
    & input_params => propagator_input_elements
  use unit_test_framework, only : unit_test_type

  implicit none

  private

  public :: propagators_test_driver

  real(dp), parameter :: tol = 1e-11_dp
  real(dp), parameter :: dt = 0.1_dp
  integer(i32), parameter :: order_Taylor = 4
  integer(i32), parameter :: dim = 2 ! Effective dimension to take into account
  integer(i32), parameter :: n_states = 2
  integer(i32), parameter :: n_kpoints = 2
  integer(i32), parameter :: n_eigvecs_houston = 2
  
  complex(dp), parameter :: H_minus_dt_hermitian(dim+1, dim+1, n_kpoints) = &
    reshape( [0.6_dp*herm(1:dim+1, 1:dim+1), 2*herm(1:dim+1, 1:dim+1)], [dim+1, dim+1, n_kpoints] )
  complex(dp), parameter :: H_0_hermitian(dim+1, dim+1, n_kpoints) = &
    reshape( [herm(1:dim+1, 1:dim+1), herm(1:dim+1, 1:dim+1)], [dim+1, dim+1, n_kpoints] )
  complex(dp), parameter :: H_dt_hermitian(dim+1, dim+1, n_kpoints) = &
    reshape( [1.4_dp*herm(1:dim+1, 1:dim+1), 0*herm(1:dim+1, 1:dim+1)], [dim+1, dim+1, n_kpoints] )
  complex(dp), parameter :: H_0_nonhermitian(dim+1, dim+1, n_kpoints) = &
    reshape( c5x7, [dim+1, dim+1, n_kpoints] )
  complex(dp), parameter :: H_dt_nonhermitian(dim+1, dim+1, n_kpoints) = &
    reshape( c7x5, [dim+1, dim+1, n_kpoints] )
  complex(dp), parameter :: S(dim+1, dim+1, n_kpoints) = &
    reshape( [pos(1:dim+1, 1:dim+1), pos(1:dim+1, 1:dim+1)], [dim+1, dim+1, n_kpoints] )
  complex(dp), parameter :: x_0(dim+1, n_states, n_kpoints) = &
    reshape( 0.1_dp*c5x5, [dim+1, n_states, n_kpoints] )

contains

  subroutine propagators_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional, intent(in) :: kill_on_failure
    
    type(unit_test_type) :: test_report

    character(len=*), parameter :: module_tested = 'propagators'

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Run and assert tests
    call test_propagator( 'SE', test_report )
    call test_propagator( 'EMR', test_report )
    call test_propagator( 'AETRS', test_report )
    call test_propagator( 'CFM4', test_report )
    call test_propagator( 'RK4', test_report )
    call test_propagator( 'EH', test_report )
    call test_propagator( 'EHM', test_report )

    ! report results
    if ( present( kill_on_failure ) ) then
      call test_report%report( module_tested, kill_on_failure )
    else
      call test_report%report( module_tested )
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine

  
  subroutine test_propagator( method, test_report )
    !> Name of the propagator to be tested
    character(len=*), intent(in) :: method
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    class(propagator), allocatable :: prop
    character(len=:), allocatable :: test_identifier
    complex(dp), allocatable :: psi(:, :, :)
    complex(dp) :: x_expected(dim+1, n_states, n_kpoints)
    integer(i32), parameter :: dims(n_kpoints) = spread( dim, dim=1, ncopies=n_kpoints )
    integer(i32) :: ik

    test_identifier = 'test_' // method //'_propagator'

    ! 1st test: Hermitian, given H_0 and H_dt
    call create_propagator( prop, input_params( method, dt, order_Taylor, tol, n_eigvecs_houston ), .true. )

    psi = x_0
    x_expected = x_0
    call prop%evolve( list_of_H_dt=H_dt_hermitian, list_of_H_0=H_0_hermitian, list_of_S=S, psi=psi, dims=dims )
    do ik = 1, n_kpoints
      call propagator_method_to_test( method, .true., H_dt_hermitian(:, :, ik), H_0_hermitian(:, :, ik), x_expected(:, :, ik) )
    end do
    call test_report%assert( all_close( psi , x_expected, tol ), test_identifier//' - 1st test failed.')

    ! 2nd test: Hermitian, given H_0 and H_minus_dt
    psi = x_0
    x_expected = x_0
    call prop%evolve( list_of_H_minus_dt=H_minus_dt_hermitian, list_of_H_0=H_0_hermitian, list_of_S=S, psi=psi, dims=dims )
    do ik = 1, n_kpoints
      call propagator_method_to_test( method, .true., H_dt_hermitian(:, :, ik), H_0_hermitian(:, :, ik), x_expected(:, :, ik) )
    end do
    call test_report%assert( all_close( psi , x_expected, tol ), test_identifier//' - 2nd test failed.')

    ! 3rd test: non-Hermitian
    ! only for methods based on Taylor expansion
    select case( trim( method ) ) 
      case( 'SE', 'EMR', 'AETRS', 'CFM4' )
        deallocate( prop )
        call create_propagator( prop, input_params( method, dt, order_Taylor, tol, n_eigvecs_houston ), .false. )

        psi = x_0
        x_expected = x_0
        call prop%evolve( list_of_H_dt=H_dt_nonhermitian, list_of_H_0=H_0_nonhermitian, list_of_S=S, psi=psi, dims=dims )
        do ik = 1, n_kpoints
          call propagator_method_to_test( method, .false., H_dt_nonhermitian(:, :, ik), H_0_nonhermitian(:, :, ik), x_expected(:, :, ik) )
        end do

        call test_report%assert( all_close( a=psi , b=x_expected, tol=tol ), &
          message=test_identifier//' - 3rd test failed.')
    end select
  end subroutine


  subroutine propagator_method_to_test( method, is_hermitian, H_dt, H_0, x )
    !> Name of the propagator to be tested
    character(len=*), intent(in) :: method
    !> If `.true.`, choose the exponential that assumes hermitian matrices
    logical, intent(in) :: is_hermitian
    !> Hamiltonian at time \(\Delta t\)
    complex(dp), contiguous, intent(in) :: H_dt(:, :)
    !> Hamiltonian at time \(0\)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    !> In: wavefunction coefficients at time \(0\)
    !> Out: wavefunction coefficients at time \(\Delta t\)
    complex(dp), contiguous, intent(inout) :: x(:, :)

    procedure(exp_general), pointer :: exp_operator 
    real(dp), parameter :: f1 =  0.5_dp - sqrt(3._dp)/6 ! 1/2 - sqrt(3)/6
    real(dp), parameter :: f2 =  0.5_dp + sqrt(3._dp)/6 ! 1/2 + sqrt(3)/6
    real(dp), parameter :: a1 =  0.25_dp - sqrt(3._dp)/6 ! 1/4 - sqrt(3)/6
    real(dp), parameter :: a2 =  0.25_dp + sqrt(3._dp)/6 ! 1/4 + sqrt(3)/6
    real(dp), parameter :: b_0 = a1*(1-f2) + a2*(1-f1)
    real(dp), parameter :: b_dt = a1*f2 + a2*f1
    real(dp), parameter :: c_0 = a1*(1-f1) + a2*(1-f2)
    real(dp), parameter :: c_dt = a1*f1 + a2*f2

    if( is_hermitian ) then
      exp_operator => exp_hermitian
    else
      exp_operator => exp_general
    end if

    select case( trim( method ) )
      case( 'SE' )
        call exp_operator( order_Taylor, -zi*dt, H_0, S(:, :, 1), x )
      case( 'EMR' )
        call exp_operator( order_Taylor, -zi*dt, 0.5_dp*(H_0+H_dt), S(:, :, 1), x )
      case( 'AETRS' )
        call exp_operator( order_Taylor, -zi*dt, 0.5_dp*H_0, S(:, :, 1), x )
        call exp_operator( order_Taylor, -zi*dt, 0.5_dp*H_dt, S(:, :, 1), x )
      case( 'CFM4' )
        call exp_operator( order_Taylor, -zi*dt, b_0*H_0 + b_dt*H_dt, S(:, :, 1), x )
        call exp_operator( order_Taylor, -zi*dt, c_0*H_0 + c_dt*H_dt, S(:, :, 1), x )
      case( 'RK4' )
        call rk4( dt, zi, H_0, 2*H_0-H_dt, S(:, :, 1), x )
      case( 'EH' )
        call exp_houston( -zi*dt, H_0(1:dim, 1:dim), S(1:dim, 1:dim, 1), x(1:dim, :), tol, n_eigvecs_houston )
      case( 'EHM' )
        call exp_houston( -zi*dt, 0.5_dp*( H_0(1:dim, 1:dim) + H_dt(1:dim, 1:dim) ), S(1:dim, 1:dim, 1), x(1:dim, :), tol, n_eigvecs_houston )
    end select
  end subroutine

end module