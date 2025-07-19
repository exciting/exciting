! Copyright 2023 EXCITING developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied. See the License for the specific language governing
! permissions and limitations under the License.

!> @file
!> Contains elements to compute a natural cubic spline, its integrals and first derivative
module idiel_cubic_spline

    use idiel_constants,   only: i32, r64, zone, zzero 

    implicit none

    type cubic_spline_t
        !> The upper limit of each spline
        real(r64), private, allocatable :: x(:) 
        !> The splines (a,b,c,d,xj)
        complex(r64), private, allocatable :: splines(:,:)
        !> The integrate of a fragment
        complex(r64), private, allocatable :: integrals(:)
    contains
        procedure, public  :: initialize=>init_cubic_spline_t, interpolate, integrate, derivative
        procedure, private :: cubic_poly_integral
        final :: clean
    end type

contains

    !> Deallocates all
    !> @param[in] this - the obj to deallocate things from
    subroutine clean(this)

        type(cubic_spline_t), intent(inout) :: this

        if (allocated(this%x))         deallocate(this%x)
        if (allocated(this%splines))   deallocate(this%splines)
        if (allocated(this%integrals)) deallocate(this%integrals)

    end subroutine clean

    !> Returns an index pointing to the first element in the array such that element < value using binary search
    !> Adapted from std c++ library routine. 
    !> @param[in] array - the array to find the lower bound of (it must be sorted in ascending order)
    !> @param[in] value - the value to compare
    !> @return answer   - the index such that array(i) < value 
    pure function lower_bound(array, value) result(answer)

        real(r64), intent(in) :: array(:)
        real(r64), intent(in) :: value 

        integer(i32) :: answer, right, mid

        answer  = 1_i32
        right   = size(array)

        do while (answer <= right)
            
            mid = (answer + right) / 2
            
            if (array(mid) == value) then
                answer = mid
                return
            else if (array(mid) < value) then
                answer = mid + 1
            else
                right = mid - 1
            end if
            
        end do

    end function lower_bound

    !> This routine initializes the splines class; it uses a natural spline algorithm
    !> that is the second derivative is supposed to be 0 at the endpoints
    !> so, the graph of the spline is a straight line outside of the interval, but still smooth. 
    !> @param[inout] this - the spline object to initialize
    !> @param[in]    x    - the variable at which the data is taken
    !> @param[in]    y    - data that need to be interpolated
    subroutine init_cubic_spline_t(this, x, y)
        class(cubic_spline_t), intent(inout) :: this
        real(r64),    intent(in)             :: x(:)
        complex(r64), intent(in)             :: y(:)

        integer(i32) :: n    ! number of splines
        integer(i32) :: i, j ! Index
        complex(r64), allocatable :: a(:), b(:), c(:), d(:) ! The splines factors d*x**3 + c*x**2 + b*x + a
        real(r64), allocatable    :: h(:) ! The step size
        real(r64), allocatable    :: l(:), mu(:)
        complex(r64), allocatable :: alpha(:), z(:)

        ! Given a set of coordinates C with size n+1 the number of splines will be n:
        n = size(y) - 1_i32
        ! Create a new array a of size n+1 and set a_i = y_i
        allocate(a,source=y)
        ! Create new arrays b and d of size n
        allocate(b(n), d(n))
        ! Create a new array h of size n and set for
        ! h(i) = x(i+1) - x; i.e. the step size
        allocate(h(n))
        do i = 1, n
            h(i) = x(i+1) - x(i) 
        end do

        ! Create a new array alpha of size n and set for
        ! alpha(i) = 3/h(i)*(a(i+1)-a(i)( - 3/h(i-1) * (a(i)-a(i-1))
        allocate(alpha(n))
        do i = 2, n
            alpha(i) = 3.0_r64 / h(i) * (a(i + 1) - a(i)) - 3.0_r64 / h(i - 1) * (a(i) - a(i - 1))
        end do

        ! Create new arrays c,l,mu and z of size n+1
        allocate(c(n + 1), l(n + 1), mu(n + 1), z(n + 1))
        ! Set l(0) = 1 and mu(0)=z(0)=0
        l(1)  = zone
        mu(1) = zzero
        z(1)  = zzero
        ! Iterate to fill l,mu,z:
        do i = 2, n
            l(i) = 2.0_r64 * (x(i + 1) - x(i - 1)) - h(i - 1) * mu(i - 1)
            mu(i) = h(i) / l(i)
            z(i) = (alpha(i) - h(i - 1) * z(i - 1)) / l(i)
        end do
        l(n+1) = zone
        z(n+1) = zzero
        c(n+1) = zzero
        
        do i = n, 1, -1
            c(i) = z(i) - mu(i) * c(i + 1)
            b(i) = (a(i + 1) - a(i)) / h(i) - h(i) * (c(i + 1) + 2 * c(i)) / 3.0_r64
            d(i) = (c(i + 1) - c(i)) / (3 * h(i))
        end do

        allocate(this%splines(5, n))
        allocate(this%integrals(n))
        allocate(this%x(n), source=x(2:))
        do i = 1, n 
            this%splines(1, i) = a(i)
            this%splines(2, i) = b(i)
            this%splines(3, i) = c(i)
            this%splines(4, i) = d(i)
            this%splines(5, i) = x(i)
            this%integrals(i)  = this%cubic_poly_integral(a(i), b(i), c(i), d(i), x(i), x(i), x(i+1))
        end do

    end subroutine init_cubic_spline_t

    !> This function computes the cubic polynonmial integral
    !> @param[in] this - the interpolator
    !> @param[in] a - a from the cubic polynomial
    !> @param[in] b - a from the cubic polynomial
    !> @param[in] c - a from the cubic polynomial
    !> @param[in] d - a from the cubic polynomial
    !> @param[in] x0 - x0 from the cubic polynomial
    !> @param[in] xa - lower limit of the integral
    !> @param[in] xb - upper limit of the integral
    !> @result    integral - the integral between a and b
    pure function cubic_poly_integral(this, a, b, c, d, x0, xa, xb) result(integral)

        class(cubic_spline_t),  intent(in) :: this
        complex(r64), intent(in)           :: a
        complex(r64), intent(in)           :: b
        complex(r64), intent(in)           :: c
        complex(r64), intent(in)           :: d
        real(r64), intent(in)              :: x0
        real(r64), intent(in)              :: xa
        real(r64), intent(in)              :: xb

        complex(r64) :: integral

        integral = a * xb - b * x0 * xb + b / 2_r64 * xb**2 - c / 3.0_r64 * (x0 - xb)**3 + d / 4.0_r64 * (xb - x0)**4 - &
                   a * xa + b * x0 * xa - b / 2_r64 * xa**2 + c / 3.0_r64 * (x0 - xa)**3 - d / 4.0_r64 * (xa - x0)**4

    end function cubic_poly_integral

    !> This function interpolates
    !> @param[in] this - the interpolator 
    !> @param[in] new_x - the x at which the data should be interpolated
    !> @result interpolated_y - the interpolated data
    function interpolate(this, new_x) result(interpolated_y)
        
        class(cubic_spline_t), target, intent(in) :: this
        real(r64), intent(in)                     :: new_x

        complex(r64) :: interpolated_y
        complex(r64), pointer :: sp(:)
        real(r64) :: dx 

        integer(i32) :: i
#ifdef DEBUG
        ! Check if within interpolation range
        if (new_x < real(this%splines(5,1)) .or. new_x > real(this%x(size(this%x)))) then
            error stop "cubic_spline_t%interpolate: x is out of the interpolator range"
        end if
#endif
        i = lower_bound(this%x, new_x)

        sp(1:5) => this%splines(1:5,i)

        dx = new_x - real(sp(5), r64)

        interpolated_y = sp(1) + dx * (sp(2) + dx * (sp(3) + sp(4)*dx)) 

        nullify(sp) 

    end function interpolate

    !> Returns the integral over the [a,b]
    !> @param[in] this - the interpolator from which to compute the integral
    !> @param[in] low  - the low limit of the integral
    !> @param[in] up   - the upper limit of the integral
    function integrate(this, low, up) result(integral)
        
        class(cubic_spline_t), intent(in) :: this
        real(r64),             intent(in) :: low 
        real(r64),             intent(in) :: up

        complex(r64) :: integral

        integer(i32) :: i, i0, i1
        real(r64)    :: x0
        complex(r64) :: a, b, c, d

#ifdef DEBUG
        ! Check that order of integral limits is correct
        if ( up < low ) then
            error stop "cubic_spline_t%integrate: the integration limits are incorrect"
        end if

        ! Check if within interpolation range
        if (low < real(this%splines(5,1)) .or. up > real(this%x(size(this%x)))) then
            error stop "cubic_spline_t%integrate: some limits are out the interpolator range"
        end if
#endif

        ! Get where the blocks are locate
        i0 = lower_bound(this%x, low)
        i1 = lower_bound(this%x, up)

        ! Init integral value
        integral = zzero

        ! Add here the middle block ranges (we simply sum from precomputed integrals of each full piece)
        integral = sum(this%integrals(i0+1:i1-1))

        ! Process starting and ending blocks 
        ! If start and end in the same block
        if (i0 == i1) then
            a  = this%splines(1,i0)
            b  = this%splines(2,i0)
            c  = this%splines(3,i0)
            d  = this%splines(4,i0)
            x0 = real(this%splines(5,i0))
            integral = integral + this%cubic_poly_integral(a, b, c, d, x0, low, up)
        else
            ! Otherwise we need to treat each of them separately
            ! Left end
            a  = this%splines(1,i0)
            b  = this%splines(2,i0)
            c  = this%splines(3,i0)
            d  = this%splines(4,i0)
            x0 = real(this%splines(5,i0))
            integral = integral + this%cubic_poly_integral(a, b, c, d, x0, low, this%x(i0))

            ! Right end
            a  = this%splines(1,i1)
            b  = this%splines(2,i1)
            c  = this%splines(3,i1)
            d  = this%splines(4,i1)
            x0 = real(this%splines(5,i1))
            integral = integral + this%cubic_poly_integral(a, b, c, d, x0, x0, up)

        end if

    end function integrate

    !> Returns the first derivative at a desired point
    !> @param[in] this - the interpolator 
    !> @param[in] x - the x at which the derivative should be interpolated
    !> @result dydx - the derivative at x
    function derivative(this, x) result(dydx)
        
        class(cubic_spline_t), target, intent(in) :: this
        real(r64), intent(in)                     :: x

        complex(r64) :: dydx
        real(r64)    :: dx
        complex(r64), pointer :: sp(:)

        integer(i32) :: i

        i = lower_bound(this%x, x)

        if (i > size(this%splines,2)) i = i - 1_i32

        sp(1:5) => this%splines(1:5,i)

        dx = x - real(sp(5), r64)

        dydx =  sp(2)  + 2.0_r64 * sp(3) * dx + 3.0_r64 * sp(4) * dx**2

        nullify(sp) 

    end function derivative

end module idiel_cubic_spline


