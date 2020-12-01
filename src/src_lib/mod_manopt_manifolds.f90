module mod_manopt_manifolds
  use mod_manopt_matrices
  implicit none

  private

  !******************************************
  ! MANIFOLDS
  !******************************************
  type, abstract, public :: manifold
    integer :: KX                          ! number of X matrices
    integer, allocatable :: DX(:,:)        ! number of rows/columns for each matrix X
    integer :: DXO(2)                      ! first and second dimension of the array X
                                           ! as allocated in the calling subroutine
    integer :: DXI(2)                      ! first and second dimension of the array X
                                           ! as used internally

    contains
      procedure, non_overridable     :: clear       ! clear memory
      procedure, non_overridable     :: dimX        ! number of unkwowns in the optimization
      procedure( dimkX_), deferred   :: dimkX       ! number of unknowns in matrix X(k)
      procedure( inner_), deferred   :: inner       ! inner product
      procedure( norm_), deferred    :: norm        ! norm
      procedure( proj_), deferred    :: proj        ! projection
      procedure( retr_), deferred    :: retr        ! retraction
      procedure( retr_), deferred    :: mexp        ! matrix exponential
      procedure( transp_), deferred  :: transp      ! parallel transport
      procedure( proj_), deferred    :: tangent     ! tangent space
      procedure( proj_), deferred    :: egrad2rgrad ! Euclidean to Riemannian gradient
  end type manifold

  type, extends( manifold), public :: euclid_manifold
    contains
      procedure :: dimkX       => euclid_dimkX
      procedure :: inner       => euclid_inner
      procedure :: norm        => euclid_norm
      procedure :: proj        => euclid_proj
      procedure :: retr        => euclid_mexp
      procedure :: mexp        => euclid_mexp
      procedure :: transp      => euclid_transp
      procedure :: tangent     => euclid_proj
      procedure :: egrad2rgrad => euclid_proj
  end type euclid_manifold
  ! constructor
  interface euclid_manifold
    module procedure :: euclid_manifold_gen
  end interface

  type, extends( manifold), public :: stiefel_manifold
    contains
      procedure :: dimkX       => stiefel_dimkX
      procedure :: inner       => stiefel_inner
      procedure :: norm        => stiefel_norm
      procedure :: proj        => stiefel_proj
      procedure :: retr        => stiefel_retr
      procedure :: mexp        => stiefel_mexp
      procedure :: transp      => stiefel_transp
      procedure :: tangent     => stiefel_proj
      procedure :: egrad2rgrad => stiefel_proj
  end type stiefel_manifold
  ! constructor
  interface stiefel_manifold
    module procedure :: stiefel_manifold_gen
  end interface

  abstract interface
    ! dimension (number of unknowns) per matrix
    function dimkX_( M, ik) result( d)
      import manifold
      class( manifold)             :: M
      integer, intent( in)         :: ik
      integer                      :: d
    end function
    ! inner product
    function inner_( M, x, y) result( r)
      import manifold, matrix
      class( manifold) :: M
      class( matrix), intent( in)  :: x    ! matrix X
      class( matrix), intent( in)  :: y    ! matrix Y
      real(8)                      :: r
    end function
    ! norm
    function norm_( M, x) result( r)
      import manifold, matrix
      class( manifold) :: M
      class( matrix), intent( in)  :: x    ! matrix M
      real(8)                      :: r
    end function
    ! projection, tangent space and
    ! Euclidean to Riemannian gradient transformation
    subroutine proj_( M, x, u, p)
      import manifold, matrix
      class( manifold) :: M
      class( matrix), intent( in)  :: x   ! matrix X
      class( matrix), intent( in)  :: u   ! projection direction U
      class( matrix), intent( out) :: p   ! projection P
    end subroutine
    ! retraction or matrix exponential
    subroutine retr_( M, x, u, t, r)
      import manifold, matrix
      class( manifold) :: M
      real(8), intent( in)         :: t   ! scaling of U
      class( matrix), intent( in)  :: x   ! matrix X
      class( matrix), intent( in)  :: u   ! matrix to retract U
      class( matrix), intent( out) :: r   ! retraction R
    end subroutine
    ! parallel transport
    subroutine transp_( M, x, y, s, u)
      import manifold, matrix
      class( manifold) :: M
      class( matrix), intent( in)  :: x   ! matrix to transport from X
      class( matrix), intent( in)  :: y   ! matrix to transport to Y
      class( matrix), intent( in)  :: s   ! matrix to transport S
      class( matrix), intent( out) :: u   ! parallel transported matrix U
    end subroutine
  end interface

  complex(8), parameter :: zzero = cmplx( 0.d0, 0.d0, 8)
  complex(8), parameter :: zone  = cmplx( 1.d0, 0.d0, 8)

  contains

    !******************************************
    ! GENERIC MANIFOLD PROCEDURES
    !******************************************
    subroutine clear( M)
      class( manifold) :: M
      if( allocated( M%DX)) deallocate( M%DX)
      return
    end subroutine clear

    function dimX( M) result( d)
      class( manifold) :: M
      integer :: d
      integer :: ik
      d = 0
      do ik = 1, M%KX
        d = d + M%dimkX( ik)
      end do
    end function dimX

    !******************************************
    ! EUCLIDEAN MANIFOLD
    !******************************************
    ! constructor
    function euclid_manifold_gen( dxo, kx, dx) result( M)
      integer, intent( in)   :: dxo(2), kx, dx(2,kx)
      type( euclid_manifold) :: M
      
      call manifold_gen( M, dxo, kx, dx)
    end function

    function euclid_dimkX( M, ik) result( d)
      class( euclid_manifold) :: M
      integer, intent( in)    :: ik
      integer                 :: d
      d = M%DX(1,ik)*M%DX(2,ik)
    end function euclid_dimkX

    function euclid_inner( M, x, y) result( r)
      class( euclid_manifold) :: M
      class( matrix), intent( in) :: x    ! matrix X
      class( matrix), intent( in) :: y    ! matrix Y
      real(8)                     :: r

      integer :: ik

      real(8), external :: ddot
      complex(8), external :: zdotc

      r = 0.d0
      if( M%DXI(1)*M%DXI(2)*M%KX < 0) return
      if( .not. same_type_as( x, y)) &
        call error( 'euclid_inner', 'Matrices have different type.')
      select type( x)
        type is( real_matrix)
          select type( y)
            type is( real_matrix)
#ifdef USEOMP
!$omp parallel default( shared) private( ik) reduction(+:r)
!$omp do
#endif
              do ik = 1, M%KX
                r = r + ddot( M%DX(1,ik)*M%DX(2,ik), x%m(1:M%DX(1,ik),1:M%DX(2,ik),ik), 1, y%m(1:M%DX(1,ik),1:M%DX(2,ik),ik), 1)
              end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
          end select
        type is( complex_matrix)
          select type( y)
            type is( complex_matrix)
#ifdef USEOMP
!$omp parallel default( shared) private( ik) reduction(+:r)
!$omp do
#endif
              do ik = 1, M%KX
                r = r + dble( zdotc( M%DX(1,ik)*M%DX(2,ik), x%m(1:M%DX(1,ik),1:M%DX(2,ik),ik), 1, y%m(1:M%DX(1,ik),1:M%DX(2,ik),ik), 1))
              end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
          end select
      end select
    end function euclid_inner

    function euclid_norm( M, x) result( r)
      class( euclid_manifold) :: M
      class( matrix), intent( in) :: x    ! matrix X
      real(8)                     :: r

      r = dsqrt( M%inner( x, x))
    end function euclid_norm

    subroutine euclid_proj( M, x, u, p)
      class( euclid_manifold) :: M
      class( matrix), intent( in)  :: x   ! matrix X
      class( matrix), intent( in)  :: u   ! projection direction U
      class( matrix), intent( out) :: p   ! projection P

      integer :: ik

      if( .not. same_type_as( x, u)) &
        call error( 'euclid_proj', 'Input matrices have different type.')
      if( .not. same_type_as( x, p)) &
        call error( 'euclid_proj', 'Output matrixes has different type.')
      select type( u)
        type is( real_matrix)
          select type( p)
            type is( real_matrix)
              p = u
          end select
        type is( complex_matrix)
          select type( p)
            type is( complex_matrix)
              p = u
          end select
      end select
      return
    end subroutine euclid_proj

    subroutine euclid_mexp( M, x, u, t, r)
      class( euclid_manifold) :: M
      real(8), intent( in)         :: t   ! scaling of U
      class( matrix), intent( in)  :: x   ! matrix X
      class( matrix), intent( in)  :: u   ! matrix to retract U
      class( matrix), intent( out) :: r   ! retraction R

      if( .not. same_type_as( x, u)) &
        call error( 'euclid_mexp', 'Input matrices have different type.')
      if( .not. same_type_as( x, r)) &
        call error( 'euclid_mexp', 'Output matrixes has different type.')
      select type( x)
        type is( real_matrix)
          select type( u)
            type is( real_matrix)
              select type( r)
                type is( real_matrix)
                  r = x + u*t
              end select
          end select
        type is( complex_matrix)
          select type( u)
            type is( complex_matrix)
              select type( r)
                type is( complex_matrix)
                  r = x + u*t
              end select
          end select
      end select
      return
    end subroutine euclid_mexp

    subroutine euclid_transp( M, x, y, s, u)
      class( euclid_manifold) :: M
      class( matrix), intent( in)  :: x   ! matrix to transport from X
      class( matrix), intent( in)  :: y   ! matrix to transport to Y
      class( matrix), intent( in)  :: s   ! matrix to transport S
      class( matrix), intent( out) :: u   ! parallel transported matrix U

      if( .not. same_type_as( x, y)) &
        call error( 'euclid_transp', 'Input matrices have different type.')
      if( .not. same_type_as( x, s)) &
        call error( 'euclid_transp', 'Input matrices have different type.')
      if( .not. same_type_as( x, u)) &
        call error( 'euclid_transp', 'Output matrixes has different type.')
      select type( s)
        type is( real_matrix)
          select type( u)
            type is( real_matrix)
              u = s
          end select
        type is( complex_matrix)
          select type( u)
            type is( complex_matrix)
              u = s
          end select
      end select
      return
    end subroutine euclid_transp

    !******************************************
    ! STIEFEL MANIFOLD
    !******************************************
    ! constructor
    function stiefel_manifold_gen( dxo, kx, dx) result( M)
      integer, intent( in)    :: dxo(2), kx, dx(2,kx)
      type( stiefel_manifold) :: M
      
      call manifold_gen( M, dxo, kx, dx)
    end function

    function stiefel_dimkX( M, ik) result( d)
      class( stiefel_manifold) :: M
      integer, intent( in)    :: ik
      integer                 :: d
      d = M%DX(1,ik)*M%DX(2,ik) - (M%DX(2,ik)*(M%DX(2,ik)+1))/2
    end function stiefel_dimkX

    function stiefel_inner( M, x, y) result( r)
      class( stiefel_manifold) :: M
      class( matrix), intent( in) :: x    ! matrix X
      class( matrix), intent( in) :: y    ! matrix Y
      real(8)                     :: r

      integer :: ik, i
      real(8) :: rk

      real(8), external :: ddot
      complex(8), external :: zdotc

      r = 0.d0
      if( M%DXI(1)*M%DXI(2)*M%KX <= 0) return
      if( .not. same_type_as( x, y)) &
        call error( 'stiefel_inner', 'Matrices have different type.')
      select type( x)
        type is( real_matrix)
          select type( y)
            type is( real_matrix)
#ifdef USEOMP
!$omp parallel default( shared) private( ik, i, rk) reduction(+:r)
!$omp do
#endif
              do ik = 1, M%KX
                rk = 0.d0
                do i = 1, M%DX(2,ik)
                  rk = rk + ddot( M%DX(1,ik), x%m(1,i,ik), 1, y%m(1,i,ik), 1)
                end do
                r = r + rk
              end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
          end select
        type is( complex_matrix)
          select type( y)
            type is( complex_matrix)
#ifdef USEOMP
!$omp parallel default( shared) private( ik, i, rk) reduction(+:r)
!$omp do
#endif
              do ik = 1, M%KX
                rk = 0.d0
                do i = 1, M%DX(2,ik)
                  rk = rk + dble( zdotc( M%DX(1,ik), x%m(1,i,ik), 1, y%m(1,i,ik), 1))
                end do
                r = r + rk
              end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
          end select
      end select
    end function stiefel_inner

    function stiefel_norm( M, x) result( r)
      class( stiefel_manifold) :: M
      class( matrix), intent( in) :: x    ! matrix X
      real(8)                     :: r

      r = dsqrt( M%inner( x, x))
    end function stiefel_norm

    subroutine stiefel_proj( M, x, u, p)
      class( stiefel_manifold) :: M
      class( matrix), intent( in)  :: x   ! matrix X
      class( matrix), intent( in)  :: u   ! projection direction U
      class( matrix), intent( out) :: p   ! projection P

      integer :: ik
      real(8), allocatable :: rauxm(:,:)
      complex(8), allocatable :: cauxm(:,:)

      if( .not. same_type_as( x, u)) &
        call error( 'stiefel_proj', 'Input matrices have different type.')
      if( .not. same_type_as( x, p)) &
        call error( 'stiefel_proj', 'Output matrixes has different type.')
      select type( u)
        type is( real_matrix)
          select type( p)
            type is( real_matrix)
              select type( x)
                type is( real_matrix)
                  p = u
#ifdef USEOMP
!$omp parallel default( shared) private( rauxm, ik)
#endif
                  allocate( rauxm(M%DXI(2),M%DXI(2)))
#ifdef USEOMP
!$omp do
#endif
                  do ik = 1, M%KX
                    if( (M%DX(1,ik) <= 0) .or. (M%DX(2,ik) <= 0)) cycle
                    call dgemm( 't', 'n', M%DX(2,ik), M%DX(2,ik), M%DX(1,ik), 1.d0, &
                           x%m(1,1,ik), M%DXI(1), &
                           u%m(1,1,ik), M%DXI(1), 0.d0, &
                           rauxm, M%DXI(2))
                    rauxm = 0.5d0*(rauxm + transpose( rauxm))
                    call dgemm( 'n', 'n', M%DX(1,ik), M%DX(2,ik), M%DX(2,ik), -1.d0, &
                           x%m(1,1,ik), M%DXI(1), &
                           rauxm, M%DXI(2), 1.d0, &
                           p%m(1,1,ik), M%DXI(1))
                  end do
#ifdef USEOMP
!$omp end do
#endif
                  deallocate( rauxm)
#ifdef USEOMP
!$omp end parallel
#endif
              end select
          end select
        type is( complex_matrix)
          select type( p)
            type is( complex_matrix)
              select type( x)
                type is( complex_matrix)
                  p = u
#ifdef USEOMP
!$omp parallel default( shared) private( cauxm, ik)
#endif
                  allocate( cauxm(M%DXI(2),M%DXI(2)))
#ifdef USEOMP
!$omp do
#endif
                  do ik = 1, M%KX
                    if( (M%DX(1,ik) <= 0) .or. (M%DX(2,ik) <= 0)) cycle
                    call zgemm( 'c', 'n', M%DX(2,ik), M%DX(2,ik), M%DX(1,ik), zone, &
                           x%m(1,1,ik), M%DXI(1), &
                           u%m(1,1,ik), M%DXI(1), zzero, &
                           cauxm, M%DXI(2))
                    cauxm = 0.5d0*(cauxm + conjg( transpose( cauxm)))
                    call zgemm( 'n', 'n', M%DX(1,ik), M%DX(2,ik), M%DX(2,ik), -zone, &
                           x%m(1,1,ik), M%DXI(1), &
                           cauxm, M%DXI(2), zone, &
                           p%m(1,1,ik), M%DXI(1))
                  end do
#ifdef USEOMP
!$omp end do
#endif
                  deallocate( cauxm)
#ifdef USEOMP
!$omp end parallel
#endif
              end select
          end select
      end select
      return
    end subroutine stiefel_proj

    subroutine stiefel_retr( M, x, u, t, r)
      use m_linalg, only: zqr, rqr
      class( stiefel_manifold) :: M
      real(8), intent( in)         :: t   ! scaling of U
      class( matrix), intent( in)  :: x   ! matrix X
      class( matrix), intent( in)  :: u   ! matrix to retract U
      class( matrix), intent( out) :: r   ! retraction R

      integer :: ik
      real(8), allocatable :: rauxm(:,:)
      complex(8), allocatable :: cauxm(:,:)

      if( .not. same_type_as( x, u)) &
        call error( 'stiefel_retr', 'Input matrices have different type.')
      if( .not. same_type_as( x, r)) &
        call error( 'stiefel_retr', 'Output matrixes has different type.')
      select type( x)
        type is( real_matrix)
          select type( u)
            type is( real_matrix)
              select type( r)
                type is( real_matrix)
                  r = x
#ifdef USEOMP
!$omp parallel default( shared) private( rauxm, ik)
#endif
                  allocate( rauxm(M%DXI(1),M%DXI(2)))
#ifdef USEOMP
!$omp do
#endif
                  do ik = 1, M%KX
                    if( (M%DX(1,ik) .le. 0) .or. (M%DX(2,ik) .le. 0)) cycle
                    if( t .ne. 0.d0) then
                      rauxm(1:M%DX(1,ik),1:M%DX(2,ik)) = r%m(1:M%DX(1,ik),1:M%DX(2,ik),ik) + cmplx( t, 0.d0, 8)*u%m(1:M%DX(1,ik),1:M%DX(2,ik),ik)
                      call rqr( rauxm(1:M%DX(1,ik),1:M%DX(2,ik)), q=r%m(1:M%DX(1,ik),1:M%DX(2,ik),ik))
                    end if
                  end do
#ifdef USEOMP
!$omp end do
#endif
                  deallocate( rauxm)
#ifdef USEOMP
!$omp end parallel
#endif
              end select
          end select
        type is( complex_matrix)
          select type( u)
            type is( complex_matrix)
              select type( r)
                type is( complex_matrix)
                  r = x
#ifdef USEOMP
!$omp parallel default( shared) private( cauxm, ik)
#endif
                  allocate( cauxm(M%DXI(1),M%DXI(2)))
#ifdef USEOMP
!$omp do
#endif
                  do ik = 1, M%KX
                    if( (M%DX(1,ik) .le. 0) .or. (M%DX(2,ik) .le. 0)) cycle
                    if( t .ne. 0.d0) then
                      cauxm(1:M%DX(1,ik),1:M%DX(2,ik)) = r%m(1:M%DX(1,ik),1:M%DX(2,ik),ik) + cmplx( t, 0.d0, 8)*u%m(1:M%DX(1,ik),1:M%DX(2,ik),ik)
                      call zqr( cauxm(1:M%DX(1,ik),1:M%DX(2,ik)), q=r%m(1:M%DX(1,ik),1:M%DX(2,ik),ik))
                    end if
                  end do
#ifdef USEOMP
!$omp end do
#endif
                  deallocate( cauxm)
#ifdef USEOMP
!$omp end parallel
#endif
              end select
          end select
      end select
      return
    end subroutine stiefel_retr

    subroutine stiefel_mexp( M, x, u, t, r)
      use m_linalg, only: zexpm, rexpm
      class( stiefel_manifold) :: M
      real(8), intent( in)         :: t   ! scaling of U
      class( matrix), intent( in)  :: x   ! matrix X
      class( matrix), intent( in)  :: u   ! matrix to retract U
      class( matrix), intent( out) :: r   ! retraction R

      integer :: ik, i
      real(8), allocatable :: rauxm1(:,:), rauxm2(:,:), rexpm1(:,:), rexpm2(:,:)
      complex(8), allocatable :: cauxm1(:,:), cauxm2(:,:), cexpm1(:,:), cexpm2(:,:)

      if( .not. same_type_as( x, u)) &
        call error( 'stiefel_mexp', 'Input matrices have different type.')
      if( .not. same_type_as( x, r)) &
        call error( 'stiefel_mexp', 'Output matrixes has different type.')
      select type( x)
        type is( real_matrix)
          select type( u)
            type is( real_matrix)
              select type( r)
                type is( real_matrix)
#ifdef USEOMP
!$omp parallel default( shared) private( rauxm1, rauxm2, rexpm1, rexpm2, ik, i)
#endif
                  allocate( rauxm1(M%DXI(2),M%DXI(2)), rexpm1(M%DXI(2),M%DXI(2)))
                  allocate( rauxm2(2*M%DXI(2),2*M%DXI(2)), rexpm2(2*M%DXI(2),2*M%DXI(2)))
#ifdef USEOMP
!$omp do
#endif
                  do ik = 1, M%KX
                    call dgemm( 't', 'n', M%DX(2,ik), M%DX(2,ik), M%DX(1,ik), t, &
                           x%m(1,1,ik), M%DXI(1), &
                           u%m(1,1,ik), M%DXI(1), 0.d0, &
                           rauxm1, M%DXI(2))
                    rauxm2 = 0.d0
                    rauxm2(              1:M%DX(2,ik),                1:M%DX(2,ik))   = rauxm1( 1:M%DX(2,ik), 1:M%DX(2,ik))
                    rauxm2( (M%DX(2,ik)+1):2*M%DX(2,ik), (M%DX(2,ik)+1):2*M%DX(2,ik)) = rauxm1( 1:M%DX(2,ik), 1:M%DX(2,ik))
                    call dgemm( 't', 'n', M%DX(2,ik), M%DX(2,ik), M%DX(1,ik), -t*t, &
                           u%m(1,1,ik), M%DXI(1), &
                           u%m(1,1,ik), M%DXI(1), 0.d0, &
                           rauxm2(1,M%DX(2,ik)+1), 2*M%DXI(2))
                    do i = 1, M%DX(2,ik)
                      rauxm2( M%DX(2,ik)+i, i) = 1.d0
                    end do
                    call rexpm( rauxm2( 1:2*M%DX(2,ik), 1:2*M%DX(2,ik)), rexpm2( 1:2*M%DX(2,ik), 1:2*M%DX(2,ik)))
                    rexpm1 = 0.d0
                    rauxm2 = 0.d0
                    call rexpm( -rauxm1( 1:M%DX(2,ik), 1:M%DX(2,ik)), rexpm1( 1:M%DX(2,ik), 1:M%DX(2,ik)))
                    call dgemm( 'n', 'n', 2*M%DX(2,ik), M%DX(2,ik), M%DX(2,ik), 1.d0, rexpm2(1,1), 2*M%DXI(2), rexpm1(1,1), M%DXI(2), 0.d0, rauxm2, 2*M%DXI(2))
                    call dgemm( 'n', 'n', M%DX(1,ik), M%DX(2,ik), M%DX(2,ik), 1.d0, x%m(1,1,ik), M%DXI(1), rauxm2(1,1), 2*M%DXI(2), 0.d0, r%m(1,1,ik), M%DXI(1))
                    call dgemm( 'n', 'n', M%DX(1,ik), M%DX(2,ik), M%DX(2,ik), t, u%m(1,1,ik), M%DXI(1), rauxm2(M%DX(2,ik)+1,1), 2*M%DXI(2), 1.d0, r%m(1,1,ik), M%DXI(1))
                  end do
#ifdef USEOMP
!$omp end do
#endif
                  deallocate( rauxm1, rauxm2, rexpm1, rexpm2)
#ifdef USEOMP
!$omp end parallel
#endif
              end select
          end select
        type is( complex_matrix)
          select type( u)
            type is( complex_matrix)
              select type( r)
                type is( complex_matrix)
#ifdef USEOMP
!$omp parallel default( shared) private( cauxm1, cauxm2, cexpm1, cexpm2, ik, i)
#endif
                  allocate( cauxm1(M%DXI(2),M%DXI(2)), cexpm1(M%DXI(2),M%DXI(2)))
                  allocate( cauxm2(2*M%DXI(2),2*M%DXI(2)), cexpm2(2*M%DXI(2),2*M%DXI(2)))
#ifdef USEOMP
!$omp do
#endif
                  do ik = 1, M%KX
                    call zgemm( 'c', 'n', M%DX(2,ik), M%DX(2,ik), M%DX(1,ik), cmplx( t, 0.d0, 8), &
                           x%m(1,1,ik), M%DXI(1), &
                           u%m(1,1,ik), M%DXI(1), zzero, &
                           cauxm1, M%DXI(2))
                    cauxm2 = zzero
                    cauxm2(              1:M%DX(2,ik),                1:M%DX(2,ik))   = cauxm1( 1:M%DX(2,ik), 1:M%DX(2,ik))
                    cauxm2( (M%DX(2,ik)+1):2*M%DX(2,ik), (M%DX(2,ik)+1):2*M%DX(2,ik)) = cauxm1( 1:M%DX(2,ik), 1:M%DX(2,ik))
                    call zgemm( 'c', 'n', M%DX(2,ik), M%DX(2,ik), M%DX(1,ik), -cmplx( t*t, 0.d0, 8), &
                           u%m(1,1,ik), M%DXI(1), &
                           u%m(1,1,ik), M%DXI(1), zzero, &
                           cauxm2(1,M%DX(2,ik)+1), 2*M%DXI(2))
                    do i = 1, M%DX(2,ik)
                      cauxm2( M%DX(2,ik)+i, i) = zone
                    end do
                    call zexpm( cauxm2( 1:2*M%DX(2,ik), 1:2*M%DX(2,ik)), cexpm2( 1:2*M%DX(2,ik), 1:2*M%DX(2,ik)))
                    cexpm1 = zzero
                    cauxm2 = zzero
                    call zexpm( -cauxm1( 1:M%DX(2,ik), 1:M%DX(2,ik)), cexpm1( 1:M%DX(2,ik), 1:M%DX(2,ik)))
                    call zgemm( 'n', 'n', 2*M%DX(2,ik), M%DX(2,ik), M%DX(2,ik), zone, cexpm2(1,1), 2*M%DXI(2), cexpm1(1,1), M%DXI(2), zzero, cauxm2, 2*M%DXI(2))
                    call zgemm( 'n', 'n', M%DX(1,ik), M%DX(2,ik), M%DX(2,ik), zone, x%m(1,1,ik), M%DXI(1), cauxm2(1,1), 2*M%DXI(2), zzero, r%m(1,1,ik), M%DXI(1))
                    call zgemm( 'n', 'n', M%DX(1,ik), M%DX(2,ik), M%DX(2,ik), cmplx( t, 0.d0, 8), u%m(1,1,ik), M%DXI(1), cauxm2(M%DX(2,ik)+1,1), 2*M%DXI(2), zone, r%m(1,1,ik), M%DXI(1))
                  end do
#ifdef USEOMP
!$omp end do
#endif
                  deallocate( cauxm1, cauxm2, cexpm1, cexpm2)
#ifdef USEOMP
!$omp end parallel
#endif
              end select
          end select
      end select
      return
    end subroutine stiefel_mexp

    subroutine stiefel_transp( M, x, y, s, u)
      class( stiefel_manifold) :: M
      class( matrix), intent( in)  :: x   ! matrix to transport from X
      class( matrix), intent( in)  :: y   ! matrix to transport to Y
      class( matrix), intent( in)  :: s   ! matrix to transport S
      class( matrix), intent( out) :: u   ! parallel transported matrix U

      if( .not. same_type_as( x, y)) &
        call error( 'stiefel_transp', 'Input matrices have different type.')
      if( .not. same_type_as( x, s)) &
        call error( 'stiefel_transp', 'Input matrices have different type.')
      if( .not. same_type_as( x, u)) &
        call error( 'stiefel_transp', 'Output matrixes has different type.')
      call stiefel_proj( M, y, s, u)
      return
    end subroutine stiefel_transp

    !******************************************
    ! MANIFOLD CONSTRUCTORS
    !******************************************
    subroutine manifold_gen( M, dxo, kx, dx)
      class( manifold), intent( out) :: M
      integer, intent( in)           :: dxo(2), kx, dx(2,kx)

      call sanityCheck( dxo, kx, dx)
      M%KX = kx
      M%DXO = dxo
      allocate( M%DX(2,KX))
      M%DX = dx
      M%DXI(1) = max( 0, maxval( M%DX(1,:))); M%DXI(2) = max( 0, maxval( M%DX(2,:)))

      contains
        subroutine sanityCheck( dxo, kx, dx)
          integer, intent( in) :: dxo(2), kx, dx(2,kx)

          if( kx < 0) &
            call error( 'genManifold', 'The number of matrices K must not be negative.')
          if( minval( dx(1,:)) .lt. 0) &
            call error( 'genManifold', 'The number of rows in X must not be negative.')
          if( minval( dx(2,:)) .lt. 0) &
            call error( 'genManifold', 'The number of columns in X must not be negative.')
          if( any( dx(1,:) - dx(2,:) .lt. 0)) &
            call error( 'genManifold', 'The number of rows in X must not be smaller than the number columns in X.')
          if( maxval( dx(1,:)) .gt. dxo(1)) &
            call error( 'genManifold', 'The number of rows in X must not be greater than the allocated dimension.')
          if( maxval( dx(2,:)) .gt. dxo(2)) &
            call error( 'genManifold', 'The number of columns in X must not be greater than the allocated dimension.')
        end subroutine sanityCheck
    end subroutine

    !******************************************
    ! UTILITY FUNCTIONS
    !******************************************
    subroutine error( proc, msg)
      character(*), intent( in) :: proc, msg
      write(*,*)
      write(*,'("Error (mod_manopt_manifolds/",a,"): ",a)') trim(proc), trim(msg)
      stop
    end subroutine error

end module mod_manopt_manifolds
