module mod_sternheimer_u
#include "asserts.fpp"
    use precision, only: dp, i32
    use math_utils, only: is_close
    use errors_warnings, only: terminate_if_true
    use ODE_solver, only: solve_radial_second_order_ODE
    use modmpi, only: mpiglobal

    implicit none 

    private
    public :: build_inhomogeneity, build_h_tilde, build_epsilon_perturbed, compute_sternheimer_u, normalize

    interface
        subroutine fderiv(m, n, x, f, g, cf)
            use precision, only: dp, i32
            
            integer(i32), intent(in)  :: m
            integer(i32), intent(in)  :: n
            real(dp), intent(in)  :: x(n)
            real(dp), intent(in)  :: f(n)
            real(dp), intent(out) :: g(n)
            real(dp), intent(out) :: cf(3, n)
        end subroutine fderiv
    end interface

contains

    !> Build the inhomogeneity of the radial Sternheimer equation. 
    pure function build_inhomogeneity(l1, l2, r_u_mix, r_u_rad, epsilon_perturbed, nr) result(inhom)

        !> l of the unperturbed radial function
        integer(i32),   intent(in)  :: l1
        !> new l-channel due to the perturbation 
        integer(i32),   intent(in)  :: l2
        !> mixed product basis function * r
        real(dp),       intent(in)  :: r_u_mix(nr)
        !> unperturbed radial function * r
        real(dp),       intent(in)  :: r_u_rad(nr)
        !> perturbed energy value
        real(dp),       intent(in)  :: epsilon_perturbed
        !> number of radial grid points in the MT
        integer(i32),   intent(in)  :: nr

        !> inhomogeneity of the radial Sternheimer equation
        real(dp) :: inhom(nr)

        inhom(1:nr) = - r_u_mix(1:nr) * r_u_rad(1:nr)

        if (l1==l2) inhom(1:nr) = inhom(1:nr) + epsilon_perturbed * r_u_rad(1:nr)
            
    end function build_inhomogeneity

    !> Compute the perturbed energy eigenvalue for the radial Sternheimer equation
    function build_epsilon_perturbed(r_u_mix, u_rad, r, nr) result(eps_perturbed)
        !> mixed product basis function * r
        real(dp),       intent(in)  :: r_u_mix(nr)
        !> unperturbed radial function
        real(dp),       intent(in)  :: u_rad(nr)
        !> radial grid
        real(dp),       intent(in)  :: r(nr)
        !> number of gridpoints in the MT
        integer(i32),   intent(in)  :: nr

        real(dp) :: eps_perturbed
        
        ! local variables
        real(dp) :: fr(nr), gr(nr), cf(3, nr)

        fr(1:nr) = u_rad(1:nr) * r_u_mix(1:nr) * u_rad(1:nr) * r(1:nr)
        ! integrate over r
        call fderiv(-1, nr, r(1:nr), fr, gr, cf)

        ! Eq.1
        eps_perturbed = gr(nr)

    end function

    !> Build the hamiltonian of the radial Sternheimer equation without the radial derivative. 
    !> Energy eigenvalue and frequency are included. 
    pure function build_h_tilde(l, epsilon, omega, nr, vr, r) result(h_tilde)
        !> new l-channel due to the perturbation 
        integer(i32),   intent(in)  :: l
        !> unperturbed energy eigenvalue
        real(dp),       intent(in)  :: epsilon
        !> frequency
        real(dp),       intent(in)  :: omega
        !> number of radial grid points in the MT
        integer(i32),   intent(in)  :: nr
        !> effective potential 
        real(dp),       intent(in)  :: vr(nr)
        !> radial grid
        real(dp),       intent(in)  :: r(nr)

        real(dp) :: h_tilde(nr)

        h_tilde(1:nr) = (l*(l+1)) / (2.0_dp * r(1:nr)**2) + vr(1:nr) - epsilon - omega

    end function

    !> Solve the radial Sternheimer equation
    subroutine compute_sternheimer_u(nr, r, h_tilde, inhomogeneity, r_u_rad)

        !> number of radial grid points in the MT
        integer(i32), intent(in) :: nr
        !> radial grid
        real(dp), intent(in) :: r(nr)
        !> h_tilde
        real(dp), intent(in) :: h_tilde(nr)
        !> inhomogeneity
        real(dp), intent(in) :: inhomogeneity(nr)
        !> radial function response times r
        real(dp), intent(out) :: r_u_rad(nr, 2)

        ! local variables
        real(dp) :: q0(nr), q1(nr)
        real(dp) :: temp
        real(dp) :: p0_initial, p1_initial, q0_initial, q1_initial
        real(dp), parameter :: atol = 1e-15_dp
        real(dp), parameter :: rtol = 1e-8_dp

        !ToDo check this again
        p1_initial = 1.0_dp 
        q1_initial = - inhomogeneity(1)

        !ToDo adjust for higher orders of energy deriv
        CALL_ASSERT( r(1) > 0.0_dp, "compute_sternheimer_u: r(1) must be > 0.")
        temp = 1.0_dp / (r(1)**2) - 2.0_dp*h_tilde(1)


        call terminate_if_true(mpiglobal, is_close(temp, 0.0_dp, rtol, atol), &
                               "Error(compute_sternheimer_u): Division by zero is not possible." )
        
        p0_initial = p1_initial/r(1) + 2.0_dp*q1_initial
        p0_initial = p0_initial / temp
        q0_initial = p1_initial*h_tilde(1) - q1_initial/r(1)
        q0_initial = q0_initial / temp
        
        call solve_radial_second_order_ODE(.false., nr, r, h_tilde, inhomogeneity, &
             p0_initial, p1_initial, q0_initial, q1_initial, r_u_rad(1:nr, 1), r_u_rad(1:nr, 2), q0, q1) 

    end subroutine compute_sternheimer_u

    !> Normalize a radial function. 
    subroutine normalize(r_u, r, nr)

        !> r * radial function u and r * first radial derivative of u
        real(dp), intent(inout)  :: r_u(nr, 2)
        !> radial grid
        real(dp), intent(in)     :: r(nr)
        !> number of radial grid points
        integer(i32), intent(in) :: nr

        !local variables
        real(dp) :: fr(nr), gr(nr), cf(3,nr)
        real(dp) :: norm_coeff

        fr(1:nr) = r_u(1:nr, 1) ** 2
        ! integrate [u(r)*r]^2 over r
        call fderiv (-1, nr, r, fr, gr, cf)
        norm_coeff = 1.0_dp / sqrt( abs( gr(nr) ) )
        r_u(1:nr, 1:2) = r_u(1:nr, 1:2) * norm_coeff

    end subroutine normalize

end module mod_sternheimer_u
