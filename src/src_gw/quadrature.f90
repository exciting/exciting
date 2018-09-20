
module mod_quadrature

    implicit none

    private

    public  :: quadrature

    interface arth
        module procedure arth_int 
        module procedure arth_rdp
    end interface arth

contains

!!----------------------------------------------------------------------
!!
!! NAME
!!  quadrature
!!
!! FUNCTION
!!  Driver routine to perform quadratures in finite domains using different techniques.
!!  The routine improves the resolution of the grid until a given accuracy is reached
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  npts=Initial number of points, only for Gauss-Legendre. At each step this number is doubled
!!  accuracy=fractional accuracy required
!!  qopt=integer flag defining the algorithm for the quadrature:
!!  ntrial=Max number of attempts
!!    1 for Trapezoidal rule, closed, O(1/N^2) 
!!    2 for Simpson based on trapezoidal,closed, O(1/N^4)
!!    3 for Midpoint rule, open, O(1/N^2)
!!    4 for midpoint rule with cancellation of leading error, open, O(1/N^4)
!!    5 for Romber integration (closed form) and extrapolation for h-->0 (order 10 is hard-coded)
!!    6 for Romber integration with midpoint rule and extrapolation for h-->0 (order 10 is hard-coded)
!!    7 for Gauss-Legendre
!!
!! OUTPUT
!!  quad=the integral
!!  ierr=0 if quadrature converged. 
!!
recursive subroutine quadrature(func, xmin, xmax, qopt, quad, ierr, ntrial, accuracy, npts)

    ! Arguments
    interface
        function func(x)
            real(8), intent(in) :: x
            real(8)             :: func
        end function func
    end interface

    real(8),    intent(in)  :: xmin, xmax
    integer(4), intent(in)  :: qopt
    real(8),    intent(out) :: quad
    integer(4), intent(out) :: ierr
    integer(4), optional, intent(in) :: ntrial
    integer(4), optional, intent(in) :: npts
    real(8),    optional, intent(in) :: accuracy
 
    ! Local
    integer(4) :: K, KM, NT, NX, NX0, it, ix
    real(8)    :: EPS, old_st, st, old_quad, dqromb
    real(8)    :: TOL
    real(8), allocatable :: h(:), s(:)
    real(8), allocatable :: wx(:), xx(:)

    ierr = 0
    TOL  = 1d-12
    EPS  = 1d-6  ; if (PRESENT(accuracy)) EPS = accuracy
    NT   = 20    ; if (PRESENT(ntrial  )) NT  = ntrial
    quad = 0.d0
                                          
    select case (qopt)

    case (1)
        ! === Trapezoidal, closed form, O(1/N^2) 
        do it = 1, NT
            call trapezoidal(func, it, xmin, xmax, quad)
            if (it > 5) then ! Avoid spurious early convergence
                if (ABS(quad-old_quad)<EPS*ABS(old_quad) .or. (ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
            end if
            old_quad = quad
        end do

    case (2)
        ! === Extended Simpson rule based on trapezoidal O(1/N^4) ===
        do it = 1, NT
            call trapezoidal(func, it, xmin, xmax, st)
            if (it == 1) then 
                quad = st
            else
                quad = (4.d0*st-old_st)/3.d0
            end if
            if (it > 5) then ! Avoid spurious early convergence
                if (ABS(quad-old_quad)<EPS*ABS(old_quad) .or. (ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
            end if
            old_quad = quad
            old_st = st
        end do

    case (3) 
        ! === Midpoint rule, open form, O(1/N^2) ===
        do it = 1, NT
            call midpoint(func, it, xmin, xmax, quad)
            if (it > 4) then ! Avoid spurious early convergence
                if (ABS(quad-old_quad)<EPS*ABS(old_quad) .or. (ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
            end if
            old_quad = quad
        end do

    case (4) 
        ! === Midpoint rule with cancellation of leading 1/N^2 term, open form, O(1/N^4) ===
        do it = 1, NT
            call midpoint(func, it, xmin, xmax, st)
            if (it == 1) then 
                quad = st
            else
                quad = (9.d0*st-old_st)/8.d0
            end if
            if (it > 4) then ! Avoid spurious early convergence
                if (ABS(quad-old_quad)<EPS*ABS(old_quad) .or. (ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
            end if
            old_quad = quad
            old_st = st
        end do

    case (5) 
        ! === Romberg Integration, closed form ===
        K=5 ; KM=K-1 ! Order 10
        allocate(h(NT+1))
        allocate(s(NT+1))
        h = 0.d0
        s = 0.d0
        h(1) = 1.d0
        do it = 1, NT
            call trapezoidal(func, it, xmin, xmax, s(it))
            !write(std_out,*) ' romberg-trap at ',ncall,it,s(it)
            if (it >= K) then 
                call polyn_interp(h(it-KM:it), s(it-KM:it), 0.d0, quad, dqromb)
                if (ABS(dqromb) < EPS*ABS(quad)) then
                    deallocate(h)
                    deallocate(s)
                    RETURN
                end if
            end if
            s(it+1) = s(it)
            h(it+1) = 0.25d0*h(it) ! Quarter makes the extrapolation a polynomial in h^2, 
        end do  ! This is required to use the Euler-Maclaurin formula 
        deallocate(h)
        deallocate(s)

    case (6) 
        ! === Romberg Integration, closed form ===
        K=5 ; KM=K-1 ! Order 10
        allocate(h(NT+1))
        allocate(s(NT+1))
        h = 0.d0 
        s = 0.d0
        h(1) = 1.d0
        do it = 1, NT
            call midpoint(func, it, xmin, xmax, s(it))
            if (it >= K) then
                call polyn_interp(h(it-KM:it), s(it-KM:it), 0.d0, quad, dqromb)
                !write(std_out,*) quad,dqromb
                if (ABS(dqromb) < EPS*ABS(quad)) then 
                    deallocate(h)
                    deallocate(s)
                    RETURN
                end if
            end if
            s(it+1) = s(it)
            h(it+1) = h(it)/9.d0 ! factor is due to step tripling in midpoint and even error series
        end do
        deallocate(h)
        deallocate(s)

    case (7) 
        ! === Gauss-Legendre ===
        NX0 = 5 ; if (PRESENT(npts)) NX0 = npts 
        NX = NX0
        do it = 1, NT
            allocate(wx(NX))
            allocate(xx(NX))
            call coeffs_gausslegint(NX, xmin, xmax, xx, wx)
            quad = 0.d0
            do ix = 1, NX 
                quad = quad+wx(ix)*func(xx(ix))
            end do
            deallocate(wx)
            deallocate(xx)
            if (it > 1) then 
                !write(std_out,*) quad
                if (ABS(quad-old_quad)<EPS*ABS(old_quad) .or. (ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
            end if
            old_quad = quad
            NX = NX + NX0
        end do

    case default 
        write(*,*)
        write(*,'(a,i3)') 'Error(mod_quadrature::quadrature) Wrong value for qopt', qopt
        write(*,*)
        stop
    end select

    write(*,*) 'Warning(mod_quadrature::quadrature):'
    write(*,'(a,i0,2(a,es14.6))') "    Results are not converged within the given accuracy. ntrial= ",NT,"; EPS= ",EPS,"; TOL= ",TOL
    ierr = -1

end subroutine quadrature

!!----------------------------------------------------------------------
!!
!! NAME
!!  trapezoidal_ (PRIVATE)
!!
!! FUNCTION
!!  Compute the n-th stage of refinement of an extended trapezoidal rule
!!  adding 2^(n-2) additional interior point in the finite range of integration
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  nn=integer defining the refinement of the mesh, each call adds 2^(n-2) additional interior points 
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  quad=the integral at the n-th stage. 
!!
!!
!! NOTES
!!  When called with nn=1, the routine returns the crudest estimate of the integral
!!  Subsequent calls with nn=2,3,... (in that sequential order) will improve the accuracy
!!  by adding 2^(n-2) additional interior points. Note that quad should not be modified between sequential calls.
!!  Subroutine is defined as recursive to allow multi-dimensional integrations
!!

recursive subroutine trapezoidal(func, nn, xmin, xmax, quad)

    ! Arguments
    interface
        function func(x)
            real(8), intent(in) :: x
            real(8)             :: func
        end function func
    end interface
    
    integer(4), intent(in)    :: nn
    real(8),    intent(in)    :: xmin, xmax
    real(8),    intent(inout) :: quad

    ! Local variables
    integer(4) :: npt, ix
    real(8)    :: space, new, yy

    select case (nn)

        case (1)
            ! === Initial crude estimate (xmax-xmin)(f1+f2)/2 ===
            quad = 0.5d0*(xmax-xmin)*(func(xmin)+func(xmax))

        case (2:)
            ! === Add npt interior points of spacing space ===
            npt=2**(nn-2) ; space=(xmax-xmin)/npt 
            ! === The new sum is combined with the old integral to give a refined integral ===
            new = 0.d0
            yy = xmin + 0.5d0*space
            do ix = 1, npt
                new = new+func(yy)
                yy = yy+space
            end do
            quad = 0.5d0*(quad+space*new) 

        case (:0)
            write(*,*)
            write(*,'(a,i3)')'Error(mod_quadrature::trapezoidal) Wrong value for nn ', nn
            write(*,*)
            stop
    
    end select

end subroutine trapezoidal

!!----------------------------------------------------------------------
!! NAME
!!  midpoint_ (PRIVATE)
!!
!! FUNCTION
!!  This routine computes the n-th stage of refinement of an extended midpoint rule.
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  nn=integer defining the refinement of the mesh, each call adds (2/3)*3n-1 additional 
!!   interior points between xmin ans xmax
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  quad=the integral at the n-th stage. 
!!
!! NOTES
!!  When called with nn=1, the routine returns as quad the crudest estimate of the integral
!!  Subsequent calls with nn=2,3,... (in that sequential order) will improve the accuracy of quad by adding
!!  (2/3)*3n-1 additional interior points. quad should not be modified between sequential calls.
!!  Subroutine is defined as recursive to allow multi-dimensional integrations
!!

recursive subroutine midpoint(func,nn,xmin,xmax,quad)

    ! Arguments
    interface
        function func(x)
            real(8), intent(in) :: x
            real(8)             :: func
        end function func
    end interface
    integer(4), intent(in)    :: nn
    real(8),    intent(in)    :: xmin,xmax
    real(8),    intent(inout) :: quad

    !   Local variables
    integer(4) :: npt, ix
    real(8)    :: space
    real(8), allocatable :: xx(:)

    select case (nn)

        case (1)
            ! === Initial crude estimate done at the middle of the interval
            quad=(xmax-xmin)*func(0.5d0*(xmin+xmax))

        case (2:)
            ! === Add npt interior points, they alternate in spacing between space and 2*space ===
            allocate(xx(2*3**(nn-2)))
            npt = 3**(nn-2) 
            space = (xmax-xmin)/(3.d0*npt)
            xx(1:2*npt-1:2) = arth(xmin+0.5d0*space,3.d0*space,npt)
            xx(2:2*npt:2) = xx(1:2*npt-1:2)+2.d0*space
            ! === The new sum is combined with the old integral to give a refined integral ===
            quad = quad/3.d0
            do ix = 1, SIZE(xx)
                quad = quad+space*func(xx(ix)) 
            end do
            deallocate(xx)

        case (:0)
            write(*,*)
            write(*,'(a,i3)')'Error(mod_quadrature::midpoint) Wrong value for nn ', nn
            write(*,*)
            stop

    end select

end subroutine midpoint

!!----------------------------------------------------------------------
!! NAME
!!  polyn_interp
!!
!! FUNCTION
!!  Given arrays xa and ya of length N, and given a value x, return a value y, and an error estimate dy. 
!!  If P(x) is the polynomial of degree N-1 such that P(xai)=yai, i=1,...,N, then the returned value y=P(x).
!!
!! INPUTS
!!  xa(:) = abscissas in ascending order
!!  ya(:) = ordinates
!!  x = the point where the set of data has to be interpolated
!!
!! OUTPUT
!!  y = the interpolated value
!!  dy = error estimate
!!
!! NOTES
!!  Based on the polint routine reported in Numerical Recipies
!!
subroutine polyn_interp(xa, ya, x, y, dy)

    ! Arguments
    real(8), intent(in)  :: xa(:), ya(:)
    real(8), intent(in)  :: x
    real(8), intent(out) :: y, dy
    ! Local variables
    integer(4) :: m, n, ns, ns_(1)
    real(8), dimension(SIZE(xa)) :: c, d, den, ho

    ! === Initialize the tables of c and d ===
    c(:) = ya(:); d(:) = ya(:); ho(:) = xa(:)-x
    
    ! === Find closest table entry and initial approximation to y ===
    ns_ = minloc(ABS(x-xa))
    ns  = ns_(1)
    y   = ya(ns)
    ns  = ns-1
    
    ! === For each column of the tableau loop over current c and d and up-date them ===
    do m = 1, n-1
        den(1:n-m) = ho(1:n-m)-ho(1+m:n)
        if (ANY(den(1:n-m)==0.d0)) then
            write(*,*)
            write(*,*) 'Error(mod_quadrature::polyn_interp) Two input xa are identical'
            write(*,*)
            stop
        end if

        den(1:n-m) = (c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m) = ho(1+m:n)*den(1:n-m) ! Update c and d 
        c(1:n-m) = ho(1:n-m)*den(1:n-m)

        if (2*ns < n-m) then ! Now decide which correction, c or d, we want to add to the 
            dy = c(ns+1)     ! accumulating value of y, The last dy added is the error indication.
        else
            dy = d(ns)
            ns = ns-1
        end if

        y=y+dy
    end do

end subroutine polyn_interp

!!----------------------------------------------------------------------
!!
!! Compute the coefficients (supports and weights)
!! for Gauss-Legendre integration.
!!
!! Input:
!! n = order of integration
!! xmin = lower bound of integration
!! xmax = upper bound of integration
!!
!! Output:
!! x(n) = array of support points
!! weights(n) = array of integration weights
!!
subroutine coeffs_gausslegint(n, xmin, xmax, x, weights)
    implicit none
    ! Arguments
    integer(4), intent(in)  :: n 
    real(8),    intent(in)  :: xmin, xmax
    real(8),    intent(out) :: x(n)
    real(8),    intent(out) :: weights(n)
    ! Local
    integer(4) :: i, j
    real(8)    :: tol, pi
    real(8)    :: xl, xmean, z, z1, p1, p2, p3, pp

    tol = 1.d-13
    pi  = 4.d0*atan(1.d0)

    xl    = (xmax-xmin)*0.5d0
    xmean = (xmax+xmin)*0.5d0

    do i = 1, (n+1)/2
        z = cos(pi*(i-0.25d0)/(n+0.5d0))
        do
            p1 = 1.d0
            p2 = 0.d0
            do j = 1, n
                p3 = p2
                p2 = p1
                p1 = ((2.d0*j - 1.d0)*z*p2 - (j-1.d0)*p3)/j
            end do
            pp = n*(p2-z*p1)/(1.0d0-z**2)
            z1 = z
            z  = z1-p1/pp
            if (abs(z-z1) < tol) exit
        end do
        x(i)           = xmean-xl*z
        x(n+1-i)       = xmean+xl*z
        weights(i)     = 2.d0*xl/((1.d0-z**2)*pp**2)
        weights(n+1-i) = weights(i)
    enddo

end subroutine

!!----------------------------------------------------------------------
!! NAME
!!  arth_int
!!
!! FUNCTION
!!  Returns an array of length nn containing an arithmetic progression whose
!!  starting value is start and whose step is step. 
!!
!! INPUTS
!!  start=initial point
!!  step=the increment
!!  nn=the number of points
!!
!! OUTPUT
!!  arth(nn)=the progression
!!
function arth_int(start, step, nn)

    ! Arguments
    integer(4), intent(in) :: start, step
    integer(4), intent(in) :: nn
    integer(4)             :: arth_int(nn)

    ! Local variables
    integer(4) :: ii

    select case (nn)

        case(1:)
            arth_int(1) = start
            do ii = 2, nn
                arth_int(ii) = arth_int(ii-1)+step
            end do

        case(0) 
            RETURN

        case(:-1)
            write(*,*)
            write(*,'(a,i4)') 'Error(mod_quadrature::arth_int) Wrong value for nn ', nn
            write(*,*)
            stop
        
    end select

end function arth_int

!!----------------------------------------------------------------------
!!
function arth_rdp(start, step, nn)
    ! Arguments
    real(8),    intent(in) :: start, step
    integer(4), intent(in) :: nn
    real(8)                :: arth_rdp(nn)
    ! Local variables
    integer :: ii

    select case (nn)

        case(1:)
            arth_rdp(1) = start
            do ii = 2, nn
                arth_rdp(ii) = arth_rdp(ii-1)+step
            end do

        case(0)
            RETURN

        case(:-1)
            write(*,*)
            write(*,'(a,i4)') 'Error(mod_quadrature::arth_rdp) Wrong value for nn ', nn
            write(*,*)
            stop

    end select

end function arth_rdp

end module