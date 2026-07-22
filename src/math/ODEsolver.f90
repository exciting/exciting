!> This module contains routines that solve a radial second order ordinary differential equation
!> using the Bulirsch-Stoer and 4th-order Runge-Kutta algorithms. 
module ODE_solver
    use precision, only: i32, dp
    use modmpi, only: terminate_if_false
    use mod_spline_coefficients, only: get_property_from_spline_coefficients

    implicit none

    private
    public :: solve_radial_second_order_ODE

    contains

    !> Solves a radial second order ordinary differential equation:
    !>   \begin{align*}
    !>      \big[ - \frac{1}{2} \frac{d^2}{dr^2} +\alpha \big] r \; u(r) = f_{\text{inhom}}(r),
    !>   \end{align*}
    !> by reformulating the second order ODE into two first order ODEs:
    !>   \begin{align*}
    !>      \big (I \big ) \; \frac{d}{dr} P(r) = 2 \; Q(r) + \frac{1}{r} \;  P(r)
    !>   \end{align*}
    !>   \begin{align*}
    !>      \big (II \big ) \; -\frac{d}{dr} Q(r) - \frac{Q(r)}{r} + \alpha \;  P(r) =  f_{\text{inhom}}(r)
    !>   \end{align*}
    !> with
    !>   \begin{align*}
    !>      P(r) = r \; u(r)
    !>   \end{align*}
    !> and
    !>   \begin{align*}
    !>      Q(r) = \frac{1}{2} \; r \; \frac{d}{dr} u(r).
    !>   \end{align*}
    !> These are then integrated using the following equation to relate the expressions at two neighboring grid points:
    !>   \begin{align*}
    !>      P(r_{i}) = P(r_{i-1}) + \bigg [ \frac{d}{dr}  P(r) \bigg |_{r_{i-1}} + 
    !>      \frac{d}{dr}  P(r)\bigg|_{r_{i}} \bigg ] \; \frac{r_{i}-r_{i-1}}{2}
    !>   \end{align*}
    !> and
    !>   \begin{align*}    
    !>      Q(r_{i}) = Q(r_{i-1}) + \bigg [ \frac{d}{dr}  Q(r) \bigg |_{r_{i-1}} + 
    !>      \frac{d}{dr}  Q(r)\bigg|_{r_{i}} \bigg ] \; \frac{r_{i}-r_{i-1}}{2}
    !>   \end{align*}
    !> and the Bulirsch-Stoer method for step size control.
    subroutine solve_radial_second_order_ODE (inwards, nr, r, alpha, inhomogeneity, p0_initial, p1_initial, q0_initial, q1_initial, p0, p1, q0, q1)  
        !> if true, integration from effective infinity inwards, else integration from 0 outwards
        logical, intent(in) :: inwards
        !> number of radial grid points
        integer(i32), intent(in) :: nr
        !> radial grid
        real(dp), intent(in) :: r(nr)
        !> alpha
        real(dp), intent(in) :: alpha(nr)
        !> inhomogeneity 
        real(dp), intent(in) :: inhomogeneity(nr)
        !> boundary condition for p0
        real(dp), intent(in) :: p0_initial
        !> boundary condition for p1
        real(dp), intent(in) :: p1_initial
        !> boundary condition for q0
        real(dp), intent(in) :: q0_initial
        !> boundary condition for q1
        real(dp), intent(in) :: q1_initial
        !> P(r_i)
        real(dp), intent(out) :: p0(nr)
        !> d/dr P(r_i)
        real(dp), intent(out) :: p1(nr)
        !> Q(r_i)
        real(dp), intent(out) :: q0(nr)
        !> d/dr Q(r_i)
        real(dp), intent(out) :: q1(nr)

        
        ! local variables
        ! radial grid points
        integer(i32) :: ir
        ! radial grid point - 1
        integer(i32) :: ir_local
        ! alpha times r on radial mesh
        real (dp) :: r_alpha(nr)
        !> cubic spline coefficients (1,2,3) and work space (4) of alpha
        real (dp) :: cf_alpha(4,nr)
        ! inhomogeneity times r on radial mesh
        real (dp) :: r_inhomogeneity(nr)
        ! cubic spline coefficients (1,2,3) and work space (4) of inhomogeneity
        real (dp) :: cf_inhomogeneity(4,nr)
        ! maximal number of Bulirsch-Stoer steps 
        integer(i32), parameter :: steps_max = 32
        ! bound on the local truncation error
        real(dp), parameter :: error_bound = 1.0e-10_dp
        ! integration interval
        integer(i32) :: ir_initial, ir_final, direction

        call terminate_if_false(all(r(1:nr) > 0.0_dp), &
        &'ERROR(solve_radial_second_order_ODE): Radial grid contains 0 or non-positive values. All elements must be greater than 0.')
       
        r_alpha(1:nr) = r(1:nr) * alpha(1:nr)
        call spline4(nr, r, 1, r_alpha, cf_alpha)
        r_inhomogeneity(1:nr) = r(1:nr) * inhomogeneity(1:nr)
        call spline4(nr, r, 1, r_inhomogeneity, cf_inhomogeneity)
       
        call set_integration_interval(inwards, nr, ir_initial, ir_final, direction)

        ! set boundary conditions
        p0(ir_initial) = p0_initial
        p1(ir_initial) = p1_initial
        q0(ir_initial) = q0_initial
        q1(ir_initial) = q1_initial

        do ir = ir_initial+direction, ir_final, direction
            ir_local = ir-direction
            call bulirsch_stoer_method(p0(ir_local), p1(ir_local), &        ! P and dP/dr at (i-1)'th grid point
                                       & q0(ir_local), q1(ir_local), &      ! Q and dQ/dr at (i-1)'th grid point
                                       & r(ir_local),  &                ! (i-1)'th grid point
                                       & r(ir), &                   ! i'th grid point
                                       & r_inhomogeneity(ir_local), &   ! f_inhom at (i-1)'th grid point
                                       & cf_inhomogeneity(:,ir_local), &! spline coefficients of f_inhom at (i-1)'th grid point
                                       & r_alpha(ir_local), &              ! alpha times r at (i-1)'th grid point
                                       & cf_alpha(:,ir_local), &           ! spline coefficients of alpha at (i-1)'th grid point
                                       & steps_max, &               ! maximum number of steps used
                                       & error_bound, &             ! maximum error 
                                       & p0(ir), p1(ir), &          ! P and Q at i'th grid point
                                       & q0(ir), q1(ir))            ! radial derivative of P and Q at i'th grid point
        end do !ir

    end subroutine solve_radial_second_order_ODE

    !> Set interval and direction of integration. If the integration is performed from the smallest r 
    !> to the largest r (outward integration), the direction is set to +1. If the integration is performed
    !> from the largest r to the smallest r (inwards integration), the direction is set to -1.
    subroutine set_integration_interval(inwards, nr, ir_initial, ir_final, direction)
        !> if true, integration from the largest r inwards, else integration from smallest r outwards
        logical, intent(in) :: inwards
        !> number of radial grid points
        integer(i32), intent(in) :: nr
        !> index of first gridpoint for integration
        integer(i32), intent(out) :: ir_initial
        !> index of last gridpoint for integration
        integer(i32), intent(out) :: ir_final
        !> direction of integration
        integer(i32), intent(out) :: direction

        if (inwards) then
            ! integration interval
            ir_initial = nr
            ir_final = 1
            
            ! direction
            direction = -1 

        else
            ! integration interval
            ir_initial = 1
            ir_final = nr
            
            ! direction
            direction = 1

        endif

    end subroutine set_integration_interval

    !> Performs the Bulirsch-Stoer algorithm (Press, William H., et al. "Numerical Recipes" (1992))
    subroutine bulirsch_stoer_method(p0_initial, p1_initial, q0_initial, q1_initial, r1, r2, r_inhomogeneity, &
            cf_inhomogeneity, r_alpha, cf_alpha, steps_max, error_bound, &
            p0_return, p1_return, q0_return, q1_return)
        !> P(r_i-1)
        real(dp), intent(in) :: p0_initial
        !> Q(r_i-1)
        real(dp), intent(in) :: q0_initial
        !> d/dr P(r_i-1)
        real(dp), intent(in) :: p1_initial
        !> d/dr Q(r_i-1)
        real(dp), intent(in) :: q1_initial
        !> first grid point r_i-1
        real(dp), intent (in) :: r1
        !> second grid point r_i
        real(dp), intent (in) :: r2
        !> inhomogeneity at r_i times r
        real(dp), intent(in) :: r_inhomogeneity
        !> cubic spline coefficients for inhomogeneity at r_i
        real (dp), intent(in)  :: cf_inhomogeneity(4)
        !> potential at r_i times r
        real(dp), intent(in) :: r_alpha
        !> cubic spline coefficients for potential at r_i
        real (dp), intent(in)  :: cf_alpha(4)
        !> maximal number of Bulirsch-Stoer steps 
        integer(i32), intent(in) :: steps_max
        !> bound on the local truncation error
        real(dp), intent(in) :: error_bound
        !> P(r_i)
        real(dp), intent(out) :: p0_return
        !> Q(r_i)
        real(dp), intent(out) :: q0_return
        !> d/dr P(r_i)
        real(dp), intent(out) :: p1_return
        !> d/dr Q(r_i)
        real(dp), intent(out) :: q1_return

        ! local variables
        integer(i32)    :: step
        integer(i32)    :: n_subgrid_points
        real(dp)        :: r1_sub
        real(dp)        :: r2_sub
        real(dp)        :: p0_estimate(steps_max)
        real(dp)        :: q0_estimate(steps_max)
        real(dp)        :: p1_estimate(steps_max)
        real(dp)        :: q1_estimate(steps_max)
        integer(i32)    :: iter
        real(dp)        :: temp_1
        real(dp)        :: alpha_2
        real(dp)        :: inhomogeneity_2
        ! local truncation error in p
        real(dp)        :: error_p
        ! local truncation error in q
        real(dp)        :: error_q
        ! max local truncation error in p and q
        real(dp)        :: error

        step = 1
        error = 1.0_dp
        p0_estimate(:) = huge(error)
        q0_estimate(:) = huge(error)

        do while ((step <= steps_max) .and. (error>error_bound))

            n_subgrid_points = 2 * step

            call get_estimates_at_next_grid_point(p0_initial, p1_initial, q0_initial, q1_initial, n_subgrid_points, &
            & r1, r2, r_inhomogeneity, cf_inhomogeneity, r_alpha, cf_alpha, p0_estimate(steps_max-step+1), &
            &p1_estimate(steps_max-step+1), q0_estimate(steps_max-step+1), q1_estimate(steps_max-step+1)) 

            ! Use Neville's algorithm (Gragg, William B., J.SIAM Series B: Numerical Analysis 2.3, 384-403 (1965)) 
            ! to extrapolate to an infinite number of subgrid points
            do iter = 1, step-1
                temp_1 = 1.0_dp / ((step/real(step-iter, kind=dp))**2 - 1.0_dp)
                p0_estimate(steps_max-step+1+iter) = p0_estimate(steps_max-step+iter) + &
                                        (p0_estimate(steps_max-step+iter) - p0_estimate(steps_max-step+1+iter)) * temp_1
                q0_estimate(steps_max-step+1+iter) = q0_estimate(steps_max-step+iter) + &
                                        (q0_estimate(steps_max-step+iter) - q0_estimate(steps_max-step+1+iter)) * temp_1
            end do

            error_p = abs( (p0_estimate(steps_max) - p0_estimate(steps_max-1) ) / p0_estimate(steps_max))
            error_q = abs( (q0_estimate(steps_max) - q0_estimate(steps_max-1) ) / q0_estimate(steps_max))

            error = max(error_p, error_q)
            step = step + 1

        end do

        p0_return = p0_estimate(steps_max)
        q0_return = q0_estimate(steps_max)

        if (error > error_bound) then 
            call runge_kutta_method(p0_initial, p1_initial, q0_initial, q1_initial, &
                r1, r2, &
                r_inhomogeneity, cf_inhomogeneity, &
                r_alpha, cf_alpha, &
                p0_return, q0_return)
        end if


        alpha_2 = get_property_from_spline_coefficients(r_alpha, cf_alpha, r1, r2) / r2
        inhomogeneity_2 =  get_property_from_spline_coefficients(r_inhomogeneity, cf_inhomogeneity, r1, r2) / r2

        p1_return = get_p1(p0_return, q0_return, r2)
        q1_return = get_q1(p0_return, q0_return, inhomogeneity_2, alpha_2, r2)

    end subroutine bulirsch_stoer_method

    !> Perform 4th order Runge-Kutta algorithm (Kutta, Wilhelm. Teubner, 1901)
    subroutine runge_kutta_method(p0_initial, p1_initial, q0_initial, q1_initial, r1, r2, r_inhomogeneity, &
        cf_inhomogeneity, r_alpha, cf_alpha, p0, q0)
        !> P(r_i-1)
        real(dp), intent(in) :: p0_initial
        !> Q(r_i-1)
        real(dp), intent(in) :: q0_initial
        !> d/dr P(r_i-1)
        real(dp), intent(in) :: p1_initial
        !> d/dr Q(r_i-1)
        real(dp), intent(in) :: q1_initial
        !> first grid point r_i-1
        real(dp), intent (in) :: r1
        !> second grid point r_i
        real(dp), intent (in) :: r2
        !> inhomogeneity at r_i times r
        real(dp), intent(in) :: r_inhomogeneity
        !> cubic spline coefficients for inhomogeneity at r_i
        real (dp), intent(in)  :: cf_inhomogeneity(4)
        !> alpha at r_i times r
        real(dp), intent(in) :: r_alpha
        !> cubic spline coefficients for alpha at r_i
        real (dp), intent(in)  :: cf_alpha(4)
        !> P(r_i)
        real(dp), intent(out) :: p0
        !> Q(r_i)
        real(dp), intent(out) :: q0

        ! local variables
        real(dp) :: p0_1, p0_2, p0_3, q0_1, q0_2, q0_3
        real(dp) :: p1_1, p1_2, p1_3, p1_4, q1_1, q1_2, q1_3, q1_4 
        real(dp) :: delta_r, inhomogeneity, alpha

        delta_r = r2-r1

        p1_1 = p1_initial
        q1_1 = q1_initial
        p0_1 = p0_initial + p1_1 * delta_r * 0.5_dp
        q0_1 = q0_initial + q1_1 * delta_r * 0.5_dp

        alpha =  get_property_from_spline_coefficients(r_alpha, cf_alpha, r1, r1+delta_r * 0.5_dp) / (r1+delta_r * 0.5_dp)
        inhomogeneity =  get_property_from_spline_coefficients(r_inhomogeneity, cf_inhomogeneity, r1, &
                                                               r1+delta_r * 0.5_dp) / (r1+delta_r * 0.5_dp)

        p1_2 = get_p1(p0_1, q0_1, r1+delta_r * 0.5_dp)
        q1_2 = get_q1(p0_1, q0_1, inhomogeneity, alpha, r1+delta_r * 0.5_dp)
        p0_2 = p0_initial + p1_2 * delta_r * 0.5_dp
        q0_2 = q0_initial + q1_2 * delta_r * 0.5_dp

        p1_3 = get_p1(p0_2, q0_2, r1+delta_r * 0.5_dp)
        q1_3 = get_q1(p0_2, q0_2, inhomogeneity, alpha, r1+delta_r * 0.5_dp)
        p0_3 = p0_initial + p1_3 * delta_r
        q0_3 = q0_initial + q1_3 * delta_r

        alpha =  get_property_from_spline_coefficients(r_alpha, cf_alpha, r1, r2) / r2
        inhomogeneity = get_property_from_spline_coefficients(r_inhomogeneity, cf_inhomogeneity, r1, r2) / r2

        p1_4 = get_p1(p0_3, q0_3, r2)
        q1_4 = get_q1(p0_3, q0_3, inhomogeneity, alpha, r2)

        p0 = p0_initial + (p1_1 + 2*p1_2 + 2*p1_3 + p1_4) *delta_r/6.0_dp
        q0 = q0_initial + (q1_1 + 2*q1_2 + 2*q1_3 + q1_4) *delta_r/6.0_dp


    end subroutine runge_kutta_method

    !> Computes the estimates of \( P(r) \), \( \frac{d}{dr} P(r) \),  \( Q(r) \), \( \frac{d}{dr} Q(r) \) 
    !> at the next grid point, from the properties given at the previous grid point, by iterating over a given number of subgrid points.  
    pure subroutine get_estimates_at_next_grid_point(p0_initial, p1_initial, q0_initial, q1_initial, n_subgrid_points, &
        & r1, r2, r_inhomogeneity, cf_inhomogeneity, r_alpha, cf_alpha, p0, p1, q0, q1)
        !> P(r_i-1)
        real(dp), intent(in) :: p0_initial
        !> Q(r_i-1)
        real(dp), intent(in) :: q0_initial
        !> d/dr P(r_i-1)
        real(dp), intent(in) :: p1_initial
        !> d/dr Q(r_i-1)
        real(dp), intent(in) :: q1_initial
        !> number of substeps
        integer(i32), intent (in) :: n_subgrid_points
        !> first grid point r_i-1
        real(dp) , intent (in) :: r1
        !> first grid point r_i
        real(dp) , intent (in) :: r2
        !> inhomogeneity at r_i times r 
        real(dp), intent(in) :: r_inhomogeneity
        !> cubic spline coefficients for inhomogeneity at r_i
        real (dp), intent(in)  :: cf_inhomogeneity(4)
        !> alpha at r_i times r
        real(dp), intent(in) :: r_alpha
        !> cubic spline coefficients for alpha at r_i
        real(dp), intent(in)  :: cf_alpha(4)
        !> P(r_i)
        real(dp), intent(out) :: p0
        !> Q(r_i)
        real(dp), intent(out) :: q0
        !> d/dr P(r_i)
        real(dp), intent(out) :: p1
        !> d/dr Q(r_i)
        real(dp), intent(out) :: q1


        ! local variables
        integer(i32) :: iter
        ! subgrid points 
        real(dp) :: r1_sub, r2_sub
        ! P(r_i)
        real(dp) :: p0_new
        ! Q(r_i)
        real(dp) :: q0_new
        ! alpha at r_j 
        real(dp) :: alpha
        ! inhomogeneity at r_j
        real(dp) :: inhomogeneity

        r2_sub = r1

        p0 = p0_initial
        q0 = q0_initial
        p1 = p1_initial
        q1 = q1_initial

        do iter = 1, n_subgrid_points

            ! r_j-1
            r1_sub = r2_sub

            !r_j
            r2_sub = compute_next_subgrid_point(r1, r2, n_subgrid_points, r1_sub)

            alpha = get_property_from_spline_coefficients(r_alpha, cf_alpha, r1, r2_sub) / r2_sub
            inhomogeneity = get_property_from_spline_coefficients(r_inhomogeneity, cf_inhomogeneity, r1, r2_sub) / r2_sub

            p0_new = get_p0(p0, q0, p1, q1, inhomogeneity, alpha, r1_sub, r2_sub)
            q0 = get_q0(p0, q0, p1, q1, inhomogeneity, alpha, r1_sub, r2_sub)

            p0 = p0_new
            p1 = get_p1(p0, q0, r2_sub)
            q1 =  get_q1(p0, q0, inhomogeneity, alpha, r2_sub)

        end do
    end subroutine get_estimates_at_next_grid_point
 
    !> Calculate  \( P \) at \( r_{j} \):
    !>   \begin{align*}
    !>      \tilde{r} = \bigg ( 1 - \frac{\Delta r}{r_j} \bigg ) \; \bigg[ 1 + \frac{2 \alpha}{\frac{1}{r_j^2} - \frac{1}{\Delta r^2}} \bigg ]
    !>   \end{align*}
    !>   \begin{align*}
    !>     P(r_j) = \frac{1}{\tilde{r}} \; \bigg \{ P(r_{j-1}) + \Delta r \bigg [ \frac{d}{dr} P(r_{j-1}) + \frac{2}{1+\frac{\Delta r}{r_j}} \big \{ Q(r_{j-1}) + \big [ \frac{d}{dr} Q(r_{j-1}) - f_{\text{inhom}}(r_j) \big ] \Delta r \big \} \bigg ] \bigg \} 
    !>   \end{align*}
    pure function get_p0(p0, q0, p1, q1, inhomogeneity,alpha, r1_sub, r2_sub) result(p0_new)
        !> P(r_j-1)
        real(dp), intent(in) :: p0
        !> Q(r_j-1)
        real(dp), intent(in) :: q0
        !> d/dr P(r_j-1)
        real(dp), intent(in) :: p1
        !> d/dr Q(r_j-1)
        real(dp), intent(in) :: q1
        !> inhomogeneity at r_i
        real(dp), intent(in) :: inhomogeneity
        !> alpha 
        real(dp), intent(in) :: alpha
        !> first subgrid point r_j-1
        real(dp), intent(in) :: r1_sub
        !> second subgrid point r_j
        real(dp), intent(in) :: r2_sub
        !> P(r_j)
        real(dp) :: p0_new

        !local variables 
        ! mean of the subgrid points r_j-1 and r_j
        real(dp) :: mean_r
        ! 1 + mean_r/r_j
        real(dp) :: const_1
        ! temporary arrays to store intermediate results
        real(dp) :: temp_1, temp_2

        mean_r = (r2_sub - r1_sub) * 0.5_dp
        const_1 = (1 + mean_r/r2_sub)     

        temp_1 = q0 + (q1 - inhomogeneity) * mean_r
        temp_1 = temp_1 * 2* mean_r
        temp_1 = (p0 + p1 * mean_r) * const_1 + temp_1

        temp_2 = alpha * 2 + 1/r2_sub**2 
        temp_2 = 1 - temp_2 * mean_r**2

        p0_new = temp_1 / temp_2
       
    end function get_p0

    !> Calculate  \( Q \) at \( r_{j} \):
    !>   \begin{align*}
    !>      \tilde{r} = 1 + \Delta r \; \bigg ( \frac{1}{r_j}+ \frac{2\alpha}{\frac{1}{r_j} + \frac{1}{\Delta r}} \bigg)
    !>   \end{align*}
    !>   \begin{align*}
    !>     Q(r_{j}) = \frac{1}{\tilde{r}} \; \bigg \{ Q(r_{j-1}) + \Delta r \bigg [ \frac{d}{dr}Q(r_{j-1}) + \frac{1}{1 - \frac{\Delta r}{r_j}} \alpha(r_{j}) \; \big \{ P(r_{j-1}) + \Delta r \; \frac{d}{dr}P(r_{j-1}) \big \} - f_{\text{inhom}}(r_j) \bigg ] \bigg \}
    !>   \end{align*}
    pure function get_q0(p0, q0, p1, q1, inhomogeneity, alpha, r1_sub, r2_sub) result(q0_new)
        !> P(r_j-1)
        real(dp), intent(in) :: p0
        !> Q(r_j-1)
        real(dp), intent(in) :: q0
        !> d/dr P(r_j-1)
        real(dp), intent(in) :: p1
        !> d/dr Q(r_j-1)
        real(dp), intent(in) :: q1
        !> inhomogeneity at r_i
        real(dp), intent(in) :: inhomogeneity
        !> alpha
        real(dp), intent(in) :: alpha
        !> first subgrid point r_j-1
        real(dp), intent(in) :: r1_sub
        !> second subgrid point r_j
        real(dp), intent(in) :: r2_sub
        !> Q(r_j)
        real(dp) :: q0_new

        !local variables 
        ! mean of the subgrid points r_j-1 and r_j
        real(dp) :: mean_r
        ! 2 / (1 + mean_r/r_i)
        real(dp) :: const_1
        ! temporary arrays 
        real(dp) :: temp_1, temp_2

        mean_r = (r2_sub - r1_sub) * 0.5_dp
        const_1 = 1 / (1 - mean_r/r2_sub)

        temp_1 = p0 + p1 * mean_r
        temp_1 = alpha * temp_1 * const_1 
        temp_1 = q1 + temp_1 - inhomogeneity
        temp_1 = q0 + temp_1 * mean_r

        temp_2 = 2.0_dp * alpha * mean_r * const_1 - 1/r2_sub 
        temp_2 = 1.0_dp - temp_2 * mean_r

        q0_new = temp_1 / temp_2

    end function get_q0

    !> Calculate \(\frac{d}{dr} P \) at \( r_{j} \):
    !>   \begin{align*}
    !>      \frac{d}{dr} P(r) = 2 \; Q(r) + \frac{1}{r} \;  P(r).
    !>   \end{align*}
    pure function get_p1(p0, q0, r2_sub) result(p1_new)
        !> P(r_j)
        real(dp), intent(in) :: p0
        !> Q(r_j)
        real(dp), intent(in) :: q0
        !> second subgrid point
        real(dp), intent(in) :: r2_sub
        !> d/dr  P(r_j)
        real(dp) :: p1_new

        p1_new = 2.0_dp * q0 + p0/r2_sub

    end function get_p1

    !> Calculate \(\frac{d}{dr} Q \) at \( r_{j} \):
    !>   \begin{align*}
    !>      \frac{d}{dr} Q(r)  = \alpha \;  P(r) - \frac{Q(r)}{r} -  f_{\text{inhom}}(r)
    !>   \end{align*}
    pure function get_q1(p0, q0, inhomogeneity, alpha, r2_sub) result(q1_new)
        !> P(r_j)
        real(dp), intent(in) :: p0
        !> Q(r_j)
        real(dp), intent(in) :: q0
        !> inhomogeneity at r_j
        real(dp), intent(in) :: inhomogeneity
        !> alpha 
        real(dp), intent(in) :: alpha
        !> second subgrid point
        real(dp), intent(in) :: r2_sub
        !> d/dr  Q(r_j)
        real(dp) :: q1_new

        q1_new = alpha * p0 - q0/r2_sub - inhomogeneity

    end function get_q1

    !> Calculates the next subgrid point \( r_j \) from the subgrid point \( r_{j-1} \) and the gridpoints \( r_i \) and \( r_{i-1} \):
    !>   \begin{align*}
    !>      r_j = r_{j-1} \; exp \bigg(\frac{log(r_i/r_{i-1})}{N_{\text{subgrid}}} \bigg).
    !>   \end{align*}
    pure function compute_next_subgrid_point(r1, r2, n_subgrid_points, initial_subgrid_point)
        !> first grid point r_i-1
        real(dp), intent (in) :: r1
        !> second grid point r_i
        real(dp), intent (in) :: r2
        !> number of substeps
        integer(i32), intent (in) :: n_subgrid_points
        !> initial subgrid point r_j-1
        real(dp), intent (in) :: initial_subgrid_point
        !> second subgrid point r_j
        real(dp) :: compute_next_subgrid_point

        ! local variables
        real(dp) :: rmult



        rmult = exp(log(r2/r1)/real(n_subgrid_points, kind=dp))

        compute_next_subgrid_point = initial_subgrid_point * rmult
        
    end function compute_next_subgrid_point

end module ODE_solver
