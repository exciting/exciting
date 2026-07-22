module ode_solver_test
  use precision, only: dp, i32
  use math_utils, only: all_close
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use ODE_solver, only: solve_radial_second_order_ODE

  implicit none

  private
  public :: radial_second_order_ode_test_driver

contains

  !> Run tests for radial second order ODE solver
  subroutine radial_second_order_ode_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> Test report object
    type(unit_test_type) :: test_report

    ! local variables 
    integer(i32) :: nr

    nr = 200
  
    ! Initialize test object
    call test_report%init(mpiglobal)

    ! Run and assert tests
    call test_outwards_no_inhomogeneity(test_report, nr)
    call test_outwards_inhomogeneity(test_report, nr)
    call test_inwards_no_inhomogeneity(test_report, nr)
    call test_inwards_inhomogeneity(test_report, nr)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('radial_second_order_ode_solver', kill_on_failure)
    else
      call test_report%report('radial_second_order_ode_solver')
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine radial_second_order_ode_test_driver

  !> Test the radial second order ODE solver using outwards integration 
  !> and a zero inhomogeneity
  subroutine test_outwards_no_inhomogeneity(test_report, nr)
    type(unit_test_type), intent(inout) ::  test_report
    integer(i32),         intent(in)    ::  nr

    integer(i32)  ::  i
    real(dp)      ::  r(nr), alpha(nr), inhom(nr), ref_p0(nr), ref_q0(nr)
    real(dp)      ::  p0_init, p1_init, q0_init, q1_init
    real(dp)      ::  k, tol1
    real(dp)      ::  p0(nr), p1(nr), q0(nr), q1(nr), gr(nr), cf(3, nr), u(nr), du(nr)
    real(dp)      :: cB
    real(dp), parameter :: tol = 1e-10_dp
    
    do i = 1, nr
      r(i) = 0.01_dp * i    
    end do
    
    inhom(1:nr) = 0.0_dp

    ! ----------- CASE 1: alpha = 0, P(r) = r ------------
    alpha(1:nr) = 0.0_dp

    p0_init = r(1)
    p1_init = 1.0_dp
    q0_init = 0.0_dp
    q1_init = 0.0_dp

    call solve_radial_second_order_ODE (.false., nr, r, alpha, inhom, p0_init, p1_init, &
         q0_init, q1_init, p0, p1, q0, q1)

    ref_q0(1:nr) = 0.0_dp
    
    call test_report%assert(all_close(p0, r, tol), &
         'test outwards/no inhomogeneity, case 1: expected p0 = r')
    call test_report%assert(all_close(q0, ref_q0, tol), &
         'test outwards/no inhomogeneity, case 1: expected q0 = 0')

    ! ----------- CASE 2: alpha = -2, P(r) = sin(k r) ------------
    k = 2.0_dp
    alpha(1:nr) = -2

    p0_init = sin(k * r(1))
    p1_init = k * cos(k * r(1))
    q0_init = ( k * r(1) * cos(k * r(1)) - sin(k * r(1)) ) / (2.0_dp * r(1))
    q1_init = -k*k * sin(k * r(1)) / 2.0_dp - q0_init / r(1)

    call solve_radial_second_order_ODE(.false., nr, r, alpha, inhom, &
         p0_init, p1_init, q0_init, q1_init, p0, p1, q0, q1)

    do i = 1, nr
      ref_p0(i) = sin(k * r(i))
      ref_q0(i) = ( k * r(i) * cos(k * r(i)) - sin(k * r(i)) ) / (2.0_dp * r(i))
    end do

    call test_report%assert(all_close(p0, ref_p0, tol), &
         'test outwards/no inhomogeneity, case 2: expected p0 = sin(k r)')
    call test_report%assert(all_close(q0, ref_q0, tol), &
         'test outwards/no inhomogeneity, case 2: expected q0 = (k r cos(k r) - sin(k r)) / (2 r)')
    
  end subroutine test_outwards_no_inhomogeneity

  !> Test the radial second order ODE solver using outwards integration 
  !> and a non-zero inhomogeneity. 
  subroutine test_outwards_inhomogeneity(test_report, nr)
    type(unit_test_type), intent(inout) ::  test_report
    integer(i32),         intent(in)    ::  nr

    integer(i32)  ::  i
    real(dp)      ::  r(nr), alpha(nr), inhom(nr), ref_p0(nr), ref_q0(nr)
    real(dp)      ::  p0_init, p1_init, q0_init, q1_init
    real(dp)      ::  p0(nr), p1(nr), q0(nr), q1(nr)

    real(dp), parameter :: tol = 1e-10_dp
    
    do i = 1, nr
      r(i) = 0.01_dp * i    
    end do

    ! ----------- CASE 1: P(r) = r^2, alpha = 0, inhom = -1 ------------
    alpha(1:nr) = 0.0_dp
    inhom(1:nr) = -1.0_dp

    p0_init = r(1)**2
    p1_init = 2.0_dp * r(1)
    q0_init = r(1) / 2.0_dp
    q1_init = 0.5_dp

    call solve_radial_second_order_ODE(.false., nr, r, alpha, inhom, &
         p0_init, p1_init, q0_init, q1_init, p0, p1, q0, q1)

    do i = 1, nr
      ref_p0(i) = r(i)**2
      ref_q0(i) = r(i) / 2.0_dp
    end do

    call test_report%assert(all_close(p0, ref_p0, tol), &
         'test outwards/inhomogeneity, case 1: expected p0 = r^2')
    call test_report%assert(all_close(q0, ref_q0, tol), &
         'test outwards/inhomogeneity, case 1: expected q0 = r/2')

  end subroutine test_outwards_inhomogeneity

  subroutine test_inwards_no_inhomogeneity(test_report, nr)
    type(unit_test_type), intent(inout) :: test_report
    integer(i32),         intent(in)    :: nr

    integer(i32) :: i, n
    real(dp) :: r(nr), alpha(nr), inhom(nr), ref_p0(nr), ref_q0(nr), dp_dr(nr), dq_dr(nr)
    real(dp) :: p0_init, p1_init, q0_init, q1_init
    real(dp) :: p0(nr), p1(nr), q0(nr), q1(nr)
    real(dp) :: k, rmin, rmax
    
    real(dp), parameter :: tol = 1e-10_dp

    do i = 1, nr
      r(i) = 0.01_dp * i
    end do
    alpha = 0.0_dp
    inhom = 0.0_dp

    ! ----------- CASE 1: alpha = 0, P(r) = r ------------
    p0_init = r(nr)
    p1_init = 1.0_dp
    q0_init = 0.0_dp
    q1_init = 0.0_dp

    call solve_radial_second_order_ODE(.true., nr, r, alpha, inhom, &
        p0_init, p1_init, q0_init, q1_init, p0, p1, q0, q1)

    do i = 1, nr
      ref_p0(i) = r(i)
      ref_q0(i) = 0.0_dp
    end do

    call test_report%assert(all_close(p0, ref_p0, tol), &
          'test inwards/no inhomogeneity, case 1: expected p0 = r')
    call test_report%assert(all_close(q0, ref_q0, tol), &
          'test inwards/no inhomogeneity, case 1: expected q0 = 0')

    ! ----------- CASE 2: alpha = -2, P(r) = sin(k r) ------------
    k = 2.0_dp
    alpha = -0.5_dp * k*k
    inhom = 0.0_dp

    p0_init = sin(k*r(nr))
    p1_init = k*cos(k*r(nr))
    q0_init = (k*r(nr)*cos(k*r(nr)) - sin(k*r(nr))) / (2.0_dp*r(nr))
    q1_init = -k*k*sin(k*r(nr))/2.0_dp - q0_init / r(nr)

    call solve_radial_second_order_ODE(.true., nr, r, alpha, inhom, &
      p0_init, p1_init, q0_init, q1_init, p0, p1, q0, q1)

    do i = 1, nr
      ref_p0(i) = sin(k*r(i))
      ref_q0(i) = (k*r(i)*cos(k*r(i)) - sin(k*r(i))) / (2.0_dp*r(i))
    end do

    call test_report%assert(all_close(p0, ref_p0, tol), &
         'test inwards/no inhomogeneity, case 2: expected p0 = sin(k r)')
    call test_report%assert(all_close(q0, ref_q0, tol), &
         'test inwards/no inhomogeneity, case 2: expected q0 = (k r cos(k r) - sin(k r)) / (2 r)')

  end subroutine test_inwards_no_inhomogeneity

  subroutine test_inwards_inhomogeneity(test_report, nr)
    type(unit_test_type), intent(inout) :: test_report
    integer(i32),         intent(in)    :: nr

    integer(i32) :: i
    real(dp) :: r(nr), alpha(nr), inhom(nr), ref_p0(nr), ref_q0(nr)
    real(dp) :: p0_init, p1_init, q0_init, q1_init
    real(dp) :: p0(nr), p1(nr), q0(nr), q1(nr)
    real(dp), parameter :: tol = 1e-10_dp

    do i = 1, nr
      r(i) = 0.01_dp * i
    end do

    ! ----------- CASE 1: P(r) = r^2, alpha = 0, inhom = -1 ------------
    alpha = 0.0_dp
    inhom = -1.0_dp

    p0_init = r(nr)**2
    p1_init = 2.0_dp * r(nr)
    q0_init = r(nr) / 2.0_dp
    q1_init = 0.5_dp

    call solve_radial_second_order_ODE(.true., nr, r, alpha, inhom, &
         p0_init, p1_init, q0_init, q1_init, p0, p1, q0, q1)

    do i = 1, nr
      ref_p0(i) = r(i)**2
      ref_q0(i) = r(i) / 2.0_dp
    end do

    call test_report%assert(all_close(p0, ref_p0, tol), &
         'test inwards/inhomogeneity, case 1: expected p0 = r^2')
    call test_report%assert(all_close(q0, ref_q0, tol), &
         'test inwards/inhomogeneity, case 1: expected q0 = r/2')

  end subroutine test_inwards_inhomogeneity

end module ode_solver_test
