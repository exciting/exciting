module mod_manopt_solvers
  use mod_manopt_matrices
  use mod_manopt_manifolds
  implicit none

  private

  ! external function interfaces
  abstract interface
    subroutine real_cost_external( x, dxo, kx, dx, f)
      integer, intent( in)    :: dxo(2), kx, dx(2,kx)
      real(8), intent( in)    :: x(dxo(1),dxo(2),*)
      real(8), intent( out)   :: f
    end subroutine
    subroutine complex_cost_external( x, dxo, kx, dx, f)
      integer, intent( in)    :: dxo(2), kx, dx(2,kx)
      complex(8), intent( in) :: x(dxo(1),dxo(2),*)
      real(8), intent( out)   :: f
    end subroutine

    subroutine real_grad_external( x, dxo, kx, dx, gx, dgo)
      integer, intent( in)     :: dxo(2), kx, dx(2,kx), dgo(2)
      real(8), intent( in)     :: x(dxo(1),dxo(2),*)
      real(8), intent( out)    :: gx(dgo(1),dgo(2),*)
    end subroutine
    subroutine complex_grad_external( x, dxo, kx, dx, gx, dgo)
      integer, intent( in)     :: dxo(2), kx, dx(2,kx), dgo(2)
      complex(8), intent( in)  :: x(dxo(1),dxo(2),*)
      complex(8), intent( out) :: gx(dgo(1),dgo(2),*)
    end subroutine

    subroutine real_costgrad_external( x, dxo, kx, dx, f, gx, dgo, funonly)
      integer, intent( in)       :: dxo(2), kx, dx(2,kx), dgo(2)
      real(8), intent( in)       :: x(dxo(1),dxo(2),*)
      real(8), intent( out)      :: f
      real(8), intent( inout)    :: gx(dgo(1),dgo(2),*)
      logical, optional, intent( in) :: funonly
    end subroutine
    subroutine complex_costgrad_external( x, dxo, kx, dx, f, gx, dgo, funonly)
      integer, intent( in)       :: dxo(2), kx, dx(2,kx), dgo(2)
      complex(8), intent( in)    :: x(dxo(1),dxo(2),*)
      real(8), intent( out)      :: f
      complex(8), intent( inout) :: gx(dgo(1),dgo(2),*)
      logical, optional, intent( in) :: funonly
    end subroutine

    subroutine real_update_external( x, dxo, kx, dx, it, change)
      integer, intent( in)       :: dxo(2), kx, dx(2,kx), it
      real(8), intent( inout)    :: x(dxo(1),dxo(2),*)
      logical, intent( out)      :: change
    end subroutine
    subroutine complex_update_external( x, dxo, kx, dx, it, change)
      integer, intent( in)       :: dxo(2), kx, dx(2,kx), it
      complex(8), intent( inout) :: x(dxo(1),dxo(2),*)
      logical, intent( out)      :: change
    end subroutine
  end interface

  public :: real_cost_external, complex_cost_external
  public :: real_grad_external, complex_grad_external
  public :: real_costgrad_external, complex_costgrad_external
  public :: real_update_external, complex_update_external

  !******************************************
  ! SOLVERS
  !******************************************
  type, abstract, public :: solver
    ! ALGORITHMICAL PARAMETERS
    real(8) ::    eps0 = 1.d-14           ! numerical zero
    real(8) :: epsgrad = 1.d-6            ! tolerance for gradient
    real(8) ::  epsfun = -1.d123          ! tolerance for cost function
    integer ::   minit = 0                ! minimum number of iterations
    integer ::   maxit = 1000             ! maximum number of iterations
    real(8) :: minstep = 1.d-2            ! minimum step length
    integer ::  stdout = 0                ! ouptut unit (0 = no output)
    real(8) ::   noise = 1.d-3            ! random step length in severe cases

    ! OPERATIONAL VARIABLES
    real(8)                     :: cCost  ! current cost function value
    real(8)                     :: nCost  ! next cost function value
    class( matrix), allocatable :: cGX    ! current X gradient
    class( matrix), allocatable :: cSX    ! current X step direction
    class( matrix), allocatable :: nX     ! next X

    ! STATISTICS
    integer :: nCostEval = 0              ! number of cost function evaluations
    integer :: nGradEval = 0              ! number of gradient evaluations

    ! EXTERNAL PROCEDURES
    procedure( real_cost_external),        nopass, pointer :: r_cost_ext => null()
    procedure( complex_cost_external),     nopass, pointer :: c_cost_ext => null()
    procedure( real_grad_external),        nopass, pointer :: r_grad_ext => null()
    procedure( complex_grad_external),     nopass, pointer :: c_grad_ext => null()
    procedure( real_costgrad_external),    nopass, pointer :: r_costgrad_ext => null()
    procedure( complex_costgrad_external), nopass, pointer :: c_costgrad_ext => null()
    procedure( real_update_external),      nopass, pointer :: r_update_ext => null()
    procedure( complex_update_external),   nopass, pointer :: c_update_ext => null()

    contains
      procedure, non_overridable :: solve => solver_core ! solve the problem
      procedure, non_overridable :: clear                ! clear memory
      procedure, non_overridable :: cost                 ! get cost function
      procedure, non_overridable :: grad                 ! get gradient
      procedure, non_overridable :: costgrad             ! get cost function and optionally gradient
      procedure, non_overridable :: update               ! perform optional user-defined update
  end type solver

  type, extends( solver) :: cg_solver
    ! ALGORITHMICAL PARAMETERS
    character(8) ::  beta = 'fr'          ! CG update (beta) type
  end type cg_solver
  ! constructor
  interface cg_solver
    module procedure :: cg_solver_gen
  end interface

  type, extends( solver) :: lbfgs_solver
    ! ALGORITHMICAL PARAMETERS
    integer ::  memlen = 12               ! BFGS memory length

    ! OPERATIONAL VARIABLES
    class( matrix), allocatable :: SHX(:) ! history of the steps X(i+1)-X(i)
    class( matrix), allocatable :: GHX(:) ! history of the gradient differences dX(i+1)-dX(i)
    real(8), allocatable        :: RH(:)  ! history of coefficients rho
  end type lbfgs_solver
  ! constructor
  interface lbfgs_solver
    module procedure :: lbfgs_solver_gen
  end interface

  ! public solver routines
  public :: cg_solver, lbfgs_solver

  contains

    !******************************************
    ! CG SOLVER
    !******************************************
    subroutine solver_cg_core( S, M, X)
      type( cg_solver), intent( inout)    :: S
      class( manifold), intent( in)       :: M
      class( matrix), intent( inout)      :: X

      ! timing
      real(8) :: t0, t1

      integer :: it, itReset
      real(8) :: gradNorm, stepSize, minStep, snorm, alpha, beta, df0
      logical :: finish, ultimatum, change
      class( matrix), allocatable :: nGX, nSX

      if( S%stdout .gt. 0) write( S%stdout, '("# CG minimization run statistics")')
      if( S%stdout .gt. 0) write( S%stdout, '("# maximum matrix dimensions: ",4x,2i5)') M%DXI
      if( S%stdout .gt. 0) write( S%stdout, '("# number of matrices:        ",4x,i5)') M%KX
      if( S%stdout .gt. 0) then
        select type( X)
          type is( real_matrix)
            write( S%stdout, '("# dimension of optimization: ",i9)') M%dimX()
          type is( complex_matrix)
            write( S%stdout, '("# dimension of optimization: ",i9)') 2*M%dimX()
        end select
      end if
      if( S%stdout .gt. 0) write( S%stdout, '("# CG update type:            ",x,8a)') adjustr( S%beta)
      if( S%stdout .gt. 0) write( S%stdout, '("#")')
      if( S%stdout .gt. 0) write( S%stdout, '("# iteration",13x,"time",17x,"cost",11x,"|gradient|",17x,"step",17x,"beta",4x,"directional slope")')

      ! initialization
      allocate( nGX, source=X)
      allocate( nSX, source=X)
      itReset = maxval( M%DX(1,:))*maxval( M%DX(2,:))
      stepSize = 1.d0; alpha = 1.d0; beta = 0.d0
      it = 0
      finish = .false.; ultimatum = .false.
      ! we first switch off the minstep and switch it on again after a few iterations
      minStep = S%minstep; S%minstep = 0.d0
      call timesec( t0)

      ! get current cost and gradient
      if( associated( S%r_update_ext) .or. associated( S%c_update_ext)) call S%update( M, X, it, change)
      if( associated( S%r_costgrad_ext) .or. associated( S%c_costgrad_ext)) then
        call S%costgrad( M, X, S%cCost, S%cGX)
      else
        call S%cost( M, X, S%cCost)
        call S%grad( M, X, S%cGX)
      end if
      gradNorm = M%norm( S%cGX)
      ! set current search direction to negative gradient
      S%cSX = S%cGX*(-1.d0)
      call timesec( t1)
      ! write 0th interation
      if( S%stdout > 0) write( S%stdout, '(2x,i9,3(x,g20.10),21x)', advance='no') it, t1-t0, S%cCost, gradNorm

      ! start minimization loop
      MINI: do while( .true.)
        ! check termination criteria
        if( alpha .lt. S%minstep+S%eps0) then
          if( .not. ultimatum) then
            beta = 0.d0 ! reset memory
            S%minstep = 0.d0
            ultimatum = .true.
          else if( (S%minstep < S%eps0) .and. (S%noise > 0.d0)) then
            ! if we are here, probably something went wrong already and even the gradient
            ! seems not to be a descend direction. Our last try is to move a small random
            ! step away and restart the minimization.
            call randomstep( M, S, X)
            gradNorm = M%norm( S%cGX)
          else
            finish = .true.
          end if
        else
          ultimatum = .false.
        end if
        if( it >= S%maxit)        finish = .true.
        if( gradNorm < S%epsgrad) finish = .true.
        if( S%cCost < S%epsfun)   finish = .true.
        if( it < S%minit)         finish = .false.
        if( finish) exit MINI

        ! get directional derivative
        df0 = M%inner( S%cGX, S%cSX)
        snorm = M%norm( S%cSX)
        if( S%stdout > 0) write( S%stdout, '(2(x,g20.10))') beta, df0
        ! perform line search
        stepSize = snorm
        call linesearch( M, S, X, df0, stepSize)
        alpha = stepSize/snorm
        if( (alpha > minStep) .and. (S%minstep < S%eps0)) S%minstep = minStep ! turn minstep on again
        if( alpha > S%eps0) then
          ! perform custom update
          if( associated( S%r_update_ext) .or. associated( S%c_update_ext)) then
            call S%update( M, S%nX, it, change)
            if( change) then
              if( associated( S%r_cost_ext) .or. associated( S%c_cost_ext)) then
                call S%cost( M, S%nX, S%nCost)
              else
                call S%costgrad( M, S%nX, S%nCost, S%cGX, funonly=.true.)
              end if
            end if
          end if
          ! get gradient at new X
          if( associated( S%r_grad_ext) .or. associated( S%c_grad_ext)) then
            call S%grad( M, S%nX, nGX)
          else
            call S%costgrad( M, S%nX, S%nCost, nGX)
          end if
        end if
        ! get new search direction
        if( mod( it+1, itReset) == 0) then
          beta = 0.d0
          nSX = nGX*(-1.d0)
        else
          call get_direction( M, S, X, nGX, beta, nSX)
        end if
        ! update X
        X = S%nX
        ! update current values
        S%cGX = nGX; S%cSX = nSX
        gradNorm = M%norm( S%cGX)
        S%cCost = S%nCost
        it = it + 1
        call timesec( t1)
        ! write i-th iteration
        if( S%stdout > 0) write( S%stdout, '(2x,i9,4(x,g20.10))', advance='no') it, t1-t0, S%cCost, gradNorm, alpha
      end do MINI
      ! write summary
      if( S%stdout > 0) write( S%stdout, *)
      if( S%stdout > 0) write( S%stdout, '("# evaluations of cost fun. : ",i9)') S%nCostEval
      if( S%stdout > 0) write( S%stdout, '("# evaluations of gradient  : ",i9)') S%nGradEval
      deallocate( nGX, nSX)

      ! copy further statistics to solver parameters
      S%epsgrad = gradNorm
      S%maxit   = it

      contains
        subroutine get_direction( M, S, X, nGX, b, P)
          class( manifold), intent( in)    :: M
          type( cg_solver), intent( inout) :: S
          class( matrix), intent( in)      :: X, nGX
          real(8), intent( out)            :: b
          class( matrix), intent( out)     :: P

          call get_beta( M, S, X, nGX, b)
          P = nGX*(-1.d0) + S%cSX*b

          return
        end subroutine

        subroutine get_beta( M, S, X, nGX, b)
          class( manifold), intent( in)    :: M
          type( cg_solver), intent( inout) :: S
          class( matrix), intent( in)      :: X, nGX
          real(8), intent( out)            :: b

          real(8) :: n, d
          class( matrix), allocatable :: oldGrad, diff

          allocate( oldGrad, source=X)
          allocate( diff, source=X)
          ! transport old search direction to new X
          call M%transp( X, S%nX, S%cSX, oldGrad)
          S%cSX = oldGrad
          ! transport old gradient to new X
          call M%transp( X, S%nX, S%cGX, oldGrad)

          select case( trim(S%beta))
            case( 'sd') ! steepest descend
              b = 0.d0
            case( 'fr') ! Fletcher-Reeves
              d = M%inner( S%cGX, S%cGX)
              n = M%inner( nGX, nGX)
              b = max( 0.d0, n/d)
            case( 'pr') ! Polak-Ribiere
              diff = nGX - oldGrad
              d = M%inner( S%cGX, S%cGX)
              n = M%inner( nGX, diff)
              b = max( 0.d0, n/d)
            case( 'hs') ! Hestenes-Stiefel
              diff = nGX - oldGrad
              d = M%inner( diff, S%cSX)
              n = M%inner( nGX, diff)
              b = max( 0.d0, n/d)
            case( 'hz') ! Hager-Zhang
              diff = nGX - oldGrad
              d = M%inner( diff, S%cSX)
              n = M%inner( diff, nGX) - 2.d0*M%inner( diff, diff)*M%inner( S%cSX, nGX)/d
              b = max( -1.d0/( M%norm( S%cSX)*min( 1.d-2, gradNorm)), n/d)
            case default
              call error( 'solver_cg_core', 'Invalid CG update type.')
          end select

          deallocate( oldGrad, diff)
          return
        end subroutine
    end subroutine solver_cg_core

    !******************************************
    ! L-BFGS SOLVER
    !******************************************
    subroutine solver_lbfgs_core( S, M, X)
      type( lbfgs_solver), intent( inout) :: S
      class( manifold), intent( in)       :: M
      class( matrix), intent( inout)      :: X

      ! timing
      real(8) :: t0, t1

      integer :: it, im, pos
      real(8) :: gradNorm, stepSize, minStep, scaleFactor, alpha, snorm, df0, dCost, iss, isg, rho, ndiff
      logical :: finish, accepted, ultimatum, change
      class( matrix), allocatable :: gx, sx

      if( S%stdout .gt. 0) write( S%stdout, '("# L-BFGS minimization run statistics")')
      if( S%stdout .gt. 0) write( S%stdout, '("# maximum matrix dimensions: ",4x,2i5)') M%DXI
      if( S%stdout .gt. 0) write( S%stdout, '("# number of matrices:        ",4x,i5)') M%KX
      if( S%stdout .gt. 0) then
        select type( X)
          type is( real_matrix)
            write( S%stdout, '("# dimension of optimization: ",i9)') M%dimX()
          type is( complex_matrix)
            write( S%stdout, '("# dimension of optimization: ",i9)') 2*M%dimX()
        end select
      end if
      if( S%stdout .gt. 0) write( S%stdout, '("# L-BFGS memory depth:       ",4x,i5)') S%memlen
      if( S%stdout .gt. 0) write( S%stdout, '("#")')
      if( S%stdout .gt. 0) write( S%stdout, '("# iteration",13x,"time",17x,"cost",11x,"|gradient|",17x,"step",4x,"directional slope")')

      ! initialization
      allocate( gx, source=X)
      allocate( sx, source=X)
      stepSize = 1.d0; alpha = 1.d0; scaleFactor = 1.d0
      accepted = .true.; ultimatum = .false.
      it = 0; pos = 0
      finish = .false.
      ! we first switch off the minstep and switch it on again after a few iterations
      minStep = S%minstep; S%minstep = 0.d0
      call timesec( t0)

      ! get current cost and gradient
      if( associated( S%r_update_ext) .or. associated( S%c_update_ext)) call S%update( M, X, it, change)
      if( associated( S%r_costgrad_ext) .or. associated( S%c_costgrad_ext)) then
        call S%costgrad( M, X, S%cCost, S%cGX)
      else
        call S%cost( M, X, S%cCost)
        call S%grad( M, X, S%cGX)
      end if
      gradNorm = M%norm( S%cGX)
      ! set current search direction to negative gradient
      S%cSX = S%cGX*(-1.d0)
      call timesec( t1)
      ! write 0th interation
      if( S%stdout > 0) write( S%stdout, '(2x,i9,3(x,g20.10),21x)', advance='no') it, t1-t0, S%cCost, gradNorm

      ! start minimization loop
      MINI: do while( .true.)
        ! check termination criteria
        if( alpha .lt. S%minstep+S%eps0) then
          if( .not. ultimatum) then
            pos = 0 ! reset memory
            S%minstep = 0.d0
            scaleFactor = 1.d0
            ultimatum = .true.
          else if( (S%minstep < S%eps0) .and. (S%noise > 0.d0)) then
            ! if we are here, probably something went wrong already and even the gradient
            ! seems not to be a descend direction. Our last try is to move a small random
            ! step away and restart the minimization.
            call randomstep( M, S, X)
            gradNorm = M%norm( S%cGX)
          else
            finish = .true.
          end if
        else
          ultimatum = .false.
        end if
        if( it >= S%maxit)        finish = .true.
        if( gradNorm < S%epsgrad) finish = .true.
        if( S%cCost < S%epsfun)   finish = .true.
        if( it < S%minit)         finish = .false.
        if( finish) exit MINI

        ! get search direction
        call get_direction( M, S, X, scaleFactor, min( pos, S%memlen), S%cSX)
        snorm = M%norm( S%cSX)
        ! get directional derivative
        df0 = M%inner( S%cGX, S%cSX)!/snorm
        if( S%stdout > 0) write( S%stdout, '(x,g20.10)') df0
        ! perform line search
        stepSize = snorm
        call linesearch( M, S, X, df0, stepSize)
        alpha = stepSize/snorm
        if( (alpha > minStep) .and. (S%minstep < S%eps0) .and. (pos > 0)) S%minstep = minStep ! turn minstep on again
        if( alpha > S%eps0) then
          ! perform custom update
          if( associated( S%r_update_ext) .or. associated( S%c_update_ext)) then
            call S%update( M, S%nX, it, change)
            if( change) then
              if( associated( S%r_cost_ext) .or. associated( S%c_cost_ext)) then
                call S%cost( M, S%nX, S%nCost)
              else
                call S%costgrad( M, S%nX, S%nCost, S%cGX, funonly=.true.)
              end if
            end if
          end if
          ! get the iteration step
          S%cSX = S%cSX*alpha
          ! transport old gradient and direction to new X
          call M%transp( X, S%nX, S%cGX, gx)
          call M%transp( X, S%nX, S%cSX, sx)
          ! get gradient at new X
          if( associated( S%r_grad_ext) .or. associated( S%c_grad_ext)) then
            call S%grad( M, S%nX, S%cGX)
          else
            call S%costgrad( M, S%nX, S%nCost, S%cGX)
          end if
          ! get step s and gradient difference
          gx = S%cGX - gx
          iss = M%norm( sx)
          sx = sx*(1.d0/iss); gx = gx*(1.d0/iss)
          isg = M%inner( sx, gx); iss = M%inner( sx, sx)
          ! update history
          if( (iss > S%eps0) .and. (isg/iss >= 1.d-4*gradNorm)) then
            accepted = .true.
            rho = 1.d0/isg
            scaleFactor = isg/M%inner( gx, gx)
            do im = max( 1, S%memlen-pos), S%memlen-1
              call M%transp( X, S%nX, S%SHX(im+1), S%SHX(im))
              call M%transp( X, S%nX, S%GHX(im+1), S%GHX(im))
              S%RH(im) = S%RH(im+1)
            end do
            if( S%memlen > 0) then
              S%SHX(S%memlen) = sx; S%GHX(S%memlen) = gx
              S%RH(S%memlen) = rho
            end if
            if( pos < S%memlen-1) pos = pos + 1
          else
            accepted = .false.
            do im = S%memlen-pos, S%memlen
              sx = S%SHX(im); call M%transp( X, S%nX, sx, S%SHX(im))
              gx = S%GHX(im); call M%transp( X, S%nX, gx, S%GHX(im))
            end do
            !scaleFactor = 1.d0
            !pos = 0
          end if
          !write(*,'(l,g20.10)',advance='no') accepted, iss
        end if
        ndiff = M%norm( S%nX - X)
        ! update X
        X = S%nX
        ! update current values
        dCost = 2.d0*(S%nCost - S%cCost)
        !write(*,'(2g20.10)') ndiff, dCost
        gradNorm = M%norm( S%cGX)
        S%cCost = S%nCost
        it = it + 1
        call timesec( t1)
        ! write i-th iteration
        if( S%stdout > 0) write( S%stdout, '(2x,i9,4(x,g20.10))', advance='no') it, t1-t0, S%cCost, gradNorm, alpha
      end do MINI
      ! write summary
      if( S%stdout > 0) write( S%stdout, *)
      if( S%stdout > 0) write( S%stdout, '("# evaluations of cost fun. : ",i9)') S%nCostEval
      if( S%stdout > 0) write( S%stdout, '("# evaluations of gradient  : ",i9)') S%nGradEval
      deallocate( sx, gx)

      ! copy further statistics to solver parameters
      S%epsgrad = gradNorm
      S%maxit   = it

      contains
        subroutine get_direction( M, S, X, f, l, P)
          class( manifold), intent( in)    :: M
          type( lbfgs_solver), intent( in) :: S
          class( matrix), intent( in)      :: X
          real(8), intent( in)             :: f
          integer, intent( in)             :: l
          class( matrix), intent( out)     :: P

          integer :: im
          real(8) :: isq(S%memlen)

          P = S%cGX
          do im = S%memlen, S%memlen-l+1, -1
            isq(im) = S%RH(im)*M%inner( S%SHX(im), P)
            P = P - S%GHX(im)*isq(im)
          end do
          P = P*f
          do im = S%memlen-l+1, S%memlen
            isq(im) = isq(im) - S%RH(im)*M%inner( S%GHX(im), P)
            P = P + S%SHX(im)*isq(im)
          end do
          P = P*(-1.d0)
        end subroutine
    end subroutine solver_lbfgs_core

    !******************************************
    ! LINESEARCH ALGORITHMS
    !******************************************
    subroutine randomstep( M, S, X)
      class( manifold), intent( in)  :: M     ! manifold object
      class( solver), intent( inout) :: S     ! solver object
      class( matrix), intent( inout) :: X     ! current point

      integer :: ik, n
      real(8) :: twopi
      integer, allocatable :: seed(:)
      real(8), allocatable :: aux(:,:,:)
      class( matrix), allocatable :: step

      call random_seed( size=n)
      allocate( seed(n))
      seed = 0
      call random_seed( put=seed)
      twopi = 8.d0*atan(1.d0)
      allocate( step, source=S%cSX)
      select type( step)
        type is( real_matrix)
          allocate( aux(M%DXI(1),M%DXI(2),1))
          do ik = 1, M%KX
            call random_number( aux(1:M%DX(1,ik),1:M%DX(2,ik),1))
            step%m(1:M%DX(1,ik),1:M%DX(2,ik),ik) = 2.d0*aux(1:M%DX(1,ik),1:M%DX(2,ik),1) - 1.d0
          end do
          deallocate( aux)
        type is( complex_matrix)
          allocate( aux(M%DXI(1),M%DXI(2),2))
          do ik = 1, M%KX
            call random_number( aux(1:M%DX(1,ik),1:M%DX(2,ik),1))
            call random_number( aux(1:M%DX(1,ik),1:M%DX(2,ik),2))
            step%m(1:M%DX(1,ik),1:M%DX(2,ik),ik) = cmplx( aux(1:M%DX(1,ik),1:M%DX(2,ik),1)*cos( twopi*aux(1:M%DX(1,ik),1:M%DX(2,ik),2)), &
                                                          aux(1:M%DX(1,ik),1:M%DX(2,ik),1)*sin( twopi*aux(1:M%DX(1,ik),1:M%DX(2,ik),2)), 8)
          end do
          deallocate( aux)
      end select
      S%cSX = step
      deallocate( step, seed)
      call M%retr( X, S%cSX, S%noise, S%nX)
      X = S%nX
      if( associated( S%r_costgrad_ext) .or. associated( S%c_costgrad_ext)) then
        call S%costgrad( M, X, S%cCost, S%cGX)
      else
        call S%cost( M, X, S%cCost)
        call S%grad( M, X, S%cGX)
      end if
      ! set current search direction to negative gradient
      S%cSX = S%cGX*(-1.d0)

      return
    end subroutine

    subroutine linesearch( M, S, X, df0, step, verbose)
      class( manifold), intent( in)  :: M     ! manifold object
      class( solver), intent( inout) :: S     ! solver object
      class( matrix), intent( in)    :: X     ! current point
      real(8), intent( in)           :: df0   ! gradient in current search direction at current point
      real(8), intent( inout)        :: step  ! the steplength
      logical, optional, intent( in) :: verbose

      real(8), parameter :: FCONTR    = 0.5d0
      real(8), parameter :: SUFFDECR  = 1.d-2
      integer, parameter :: MAXSTEPS  = 50
      logical, parameter :: BACKTRACK = .true.
      logical, parameter :: FORCEDECR = .false.

      integer :: it
      real(8) :: decr, alpha, df, pnorm

      pnorm = M%norm( S%cSX)
      alpha = step/pnorm; it = 1
      call M%retr( X, S%cSX, alpha, S%nX)
      if( associated( S%r_cost_ext) .or. associated( S%c_cost_ext)) then
        call S%cost( M, S%nX, S%nCost)
      else
        call S%costgrad( M, S%nX, S%nCost, S%cGX, funonly=.true.)
      end if
      decr = min( 0.5d0*S%cCost, -SUFFDECR*alpha*df0)
      if( present( verbose) .and. S%stdout > 0) write(S%stdout,'(i5,2g16.6,3g20.10)') it, alpha, df0, S%cCost, S%nCost

      do while( BACKTRACK .and. ( S%nCost > S%cCost - decr))
        alpha = max( S%minstep, alpha*FCONTR)
        call M%retr( X, S%cSX, alpha, S%nX)
        if( associated( S%r_cost_ext) .or. associated( S%c_cost_ext)) then
          call S%cost( M, S%nX, S%nCost)
        else
          call S%costgrad( M, S%nX, S%nCost, S%cGX, funonly=.true.)
        end if
        it = it + 1
        decr = min( 0.5d0*S%cCost, -SUFFDECR*alpha*df0)
        if( present( verbose) .and. S%stdout > 0) write(S%stdout,'(i5,2g16.6,3g20.10)') it, alpha, df0, S%cCost, S%nCost
        if( (it > MAXSTEPS) .or. (alpha < S%minstep+S%eps0)) exit
      end do

      if( (FORCEDECR .and. (S%nCost > S%cCost)) .or. isnan( S%nCost)) then
        alpha = 0.d0
        S%nX = X
        S%nCost = S%cCost
      end if

      if( present( verbose) .and. S%stdout > 0) write(S%stdout,'(i5,2g16.6,3g20.10)') 0, 0.d0, df0, S%cCost, S%cCost

      step = alpha*pnorm
      return
    end subroutine

    !******************************************
    ! SOLVER CONSTRUCTORS
    !******************************************

    function cg_solver_gen( M, X) result( S)
      class( manifold), intent( in) :: M
      class( matrix), intent( in)   :: X
      type( cg_solver)              :: S

      select type( X)
        type is( real_matrix)
          allocate( real_matrix :: S%cGX, S%cSX, S%nX)
          S%cGX = real_matrix( X%m)
          S%cSX = real_matrix( X%m)
          S%nX  = real_matrix( X%m)
        type is( complex_matrix)
          allocate( complex_matrix :: S%cGX, S%cSX, S%nX)
          S%cGX = complex_matrix( X%m)
          S%cSX = complex_matrix( X%m)
          S%nX  = complex_matrix( X%m)
      end select
    end function

    function lbfgs_solver_gen( M, X, memlen) result( S)
      class( manifold), intent( in)  :: M
      class( matrix), intent( in)    :: X
      integer, optional, intent( in) :: memlen
      type( lbfgs_solver)            :: S

      integer :: i

      if( present( memlen)) S%memlen = max( 0, memlen)
      allocate( S%RH( S%memlen))
      S%RH = 0.d0
      select type( X)
        type is( real_matrix)
          allocate( real_matrix :: S%cGX, S%cSX, S%nX)
          S%cGX = real_matrix( X%m)*0.d0
          S%cSX = real_matrix( X%m)*0.d0
          S%nX  = real_matrix( X%m)*0.d0
          allocate( real_matrix :: S%SHX( S%memlen), S%GHX( S%memlen))
          do i = 1, S%memlen
            S%SHX(i) = real_matrix( X%m)*0.d0
            S%GHX(i) = real_matrix( X%m)*0.d0
          end do
        type is( complex_matrix)
          allocate( complex_matrix :: S%cGX, S%cSX, S%nX)
          S%cGX = complex_matrix( X%m)*0.d0
          S%cSX = complex_matrix( X%m)*0.d0
          S%nX  = complex_matrix( X%m)*0.d0
          allocate( complex_matrix :: S%SHX( S%memlen), S%GHX( S%memlen))
          do i = 1, S%memlen
            S%SHX(i) = complex_matrix( X%m)*0.d0
            S%GHX(i) = complex_matrix( X%m)*0.d0
          end do
      end select
    end function

    !******************************************
    ! UTILITY FUNCTIONS
    !******************************************
    subroutine solver_core( S, M, X)
      class( solver)                 :: S
      class( manifold), intent( in)  :: M
      class( matrix), intent( inout) :: X

      select type( S)
        type is( cg_solver)
          call solver_cg_core( S, M, X)
        type is( lbfgs_solver)
          call solver_lbfgs_core( S, M, X)
      end select
    end subroutine

    subroutine clear( S)
      class( solver) :: S
      integer :: i
      select type( S)
        type is( cg_solver)
          call S%cGX%clear()
          call S%cSX%clear()
          call S%nX%clear()
        type is( lbfgs_solver)
          call S%cGX%clear()
          call S%cSX%clear()
          call S%nX%clear()
          do i = 1, S%memlen
            call S%SHX(i)%clear()
            call S%GHX(i)%clear()
          end do
          if( allocated( S%SHX)) deallocate( S%SHX)
          if( allocated( S%GHX)) deallocate( S%GHX)
          if( allocated( S%RH)) deallocate( S%RH)
      end select
      return
    end subroutine clear

    subroutine error( proc, msg)
      character(*), intent( in) :: proc, msg
      write(*,*)
      write(*,'("Error (mod_manopt_solvers/",a,"): ",a)') trim(proc), trim(msg)
      stop
    end subroutine error

    subroutine cost( S, M, X, f)
      class( solver)                :: S
      class( manifold), intent( in) :: M
      class( matrix), intent( in)   :: X
      real(8), intent( out)         :: f

      select type( X)
        type is( real_matrix)
          if( .not. associated( S%r_cost_ext)) &
            call error( 'cost', 'For a real input matrix, the solver variable r_cost_ext must point to an external cost subroutine.')
          call S%r_cost_ext( X%m, M%DXI, M%KX, M%DX, f)
          S%nCostEval = S%nCostEval + 1
        type is( complex_matrix)
          if( .not. associated( S%c_cost_ext)) &
            call error( 'cost', 'For a complex input matrix, the solver variable c_cost_ext must point to an external cost subroutine.')
          call S%c_cost_ext( X%m, M%DXI, M%KX, M%DX, f)
          S%nCostEval = S%nCostEval + 1
      end select
      return
    end subroutine cost

    subroutine grad( S, M, X, G)
      class( solver)                :: S
      class( manifold), intent( in) :: M
      class( matrix), intent( in)   :: X
      class( matrix), intent( out)  :: G

      class( real_matrix), allocatable :: rtmp
      class( complex_matrix), allocatable :: ctmp

      if( .not. same_type_as( X, G)) &
        call error( 'grad', 'Input matrix X and gradient G must have the same type.')
      select type( X)
        type is( real_matrix)
          if( .not. associated( S%r_grad_ext)) &
            call error( 'cost', 'For a real input matrix, the solver variable r_grad_ext must point to an external gradient subroutine.')
          select type( G)
            type is( real_matrix)
              allocate( rtmp, source=x)
              call S%r_grad_ext( X%m, M%DXI, M%KX, M%DX, rtmp%m, M%DXI)
              call M%egrad2rgrad( X, rtmp, G)
              S%nGradEval = S%nGradEval + 1
              deallocate( rtmp)
          end select
        type is( complex_matrix)
          if( .not. associated( S%c_cost_ext)) &
            call error( 'cost', 'For a complex input matrix, the solver variable c_grad_ext must point to an external gradient subroutine.')
          select type( G)
            type is( complex_matrix)
              allocate( ctmp, source=x)
              call S%c_grad_ext( X%m, M%DXI, M%KX, M%DX, ctmp%m, M%DXI)
              call M%egrad2rgrad( X, ctmp, G)
              S%nGradEval = S%nGradEval + 1
              deallocate( ctmp)
          end select
      end select
      return
    end subroutine grad

    subroutine costgrad( S, M, X, f, G, funonly)
      class( solver)                 :: S
      class( manifold), intent( in)  :: M
      class( matrix), intent( in)    :: X
      real(8), intent( out)          :: f
      class( matrix), intent( inout) :: G
      logical, optional, intent( in) :: funonly

      logical :: nograd
      class( real_matrix), allocatable :: rtmp
      class( complex_matrix), allocatable :: ctmp

      nograd = .false.
      if( present( funonly)) nograd = funonly

      if( .not. same_type_as( X, G)) &
        call error( 'costgrad', 'Input matrix X and gradient G must have the same type.')
      select type( X)
        type is( real_matrix)
          if( .not. associated( S%r_costgrad_ext)) &
            call error( 'costgrad', 'For a real input matrix, the solver variable r_costgrad_ext must point to an external cost and gradient subroutine.')
          select type( G)
            type is( real_matrix)
              allocate( rtmp, source=x)
              call S%r_costgrad_ext( X%m, M%DXI, M%KX, M%DX, f, rtmp%m, M%DXI, funonly=nograd)
              S%nCostEval = S%nCostEval + 1
              if( .not. nograd) then
                call M%egrad2rgrad( X, rtmp, G)
                S%nGradEval = S%nGradEval + 1
              end if
              deallocate( rtmp)
          end select
        type is( complex_matrix)
          if( .not. associated( S%c_costgrad_ext)) &
            call error( 'costgrad', 'For a complex input matrix, the solver variable c_costgrad_ext must point to an external cost and gradient subroutine.')
          select type( G)
            type is( complex_matrix)
              allocate( ctmp, source=x)
              call S%c_costgrad_ext( X%m, M%DXI, M%KX, M%DX, f, ctmp%m, M%DXI, funonly=nograd)
              S%nCostEval = S%nCostEval + 1
              if( .not. nograd) then
                call M%egrad2rgrad( X, ctmp, G)
                S%nGradEval = S%nGradEval + 1
              end if
              deallocate( ctmp)
          end select
      end select
      return
    end subroutine costgrad

    subroutine update( S, M, X, it, change)
      class( solver)                 :: S
      class( manifold), intent( in)  :: M
      class( matrix), intent( inout) :: X
      integer, intent( in)           :: it
      logical, intent( out)          :: change

      class( real_matrix), allocatable :: rtmp
      class( complex_matrix), allocatable :: ctmp

      select type( X)
        type is( real_matrix)
          if( .not. associated( S%r_update_ext)) &
            call error( 'update', 'For a real input matrix, the solver variable r_update_ext must point to an external update subroutine.')
          allocate( rtmp, source=X)
          call S%r_update_ext( rtmp%m, M%DXI, M%KX, M%DX, it, change)
          X%m = rtmp%m
          deallocate( rtmp)
        type is( complex_matrix)
          if( .not. associated( S%c_update_ext)) &
            call error( 'update', 'For a complex input matrix, the solver variable c_update_ext must point to an external update subroutine.')
          allocate( ctmp, source=X)
          call S%c_update_ext( ctmp%m, M%DXI, M%KX, M%DX, it, change)
          X%m = ctmp%m
          deallocate( ctmp)
      end select
      return
    end subroutine update

end module mod_manopt_solvers
