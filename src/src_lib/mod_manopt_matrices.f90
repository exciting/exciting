module mod_manopt_matrices
  implicit none

  private

  !******************************************
  ! MATRICES
  !******************************************
  type, abstract, public :: matrix
    contains
      private
      procedure, public, non_overridable  :: clear
      procedure( plot_), public, deferred :: plot
      procedure( mmsub_), deferred        :: ass
      generic, public                     :: assignment(=) => ass
      procedure( mmfun_), deferred        :: add
      generic, public                     :: operator(+) => add
      procedure( mmfun_), deferred        :: sub
      generic, public                     :: operator(-) => sub
      procedure( mmfun_), deferred        :: mmmul
      procedure( mrfun_), deferred        :: mrmul
      procedure( mcfun_), deferred        :: mcmul
      generic, public                     :: operator(*) => mmmul, mrmul, mcmul
  end type matrix

  abstract interface
    subroutine plot_( x)
      import :: matrix
      class( matrix), intent( in) :: x
    end subroutine
    function mmfun_( x, y) result( z)
      import :: matrix
      class( matrix), intent( in) :: x, y
      class( matrix), allocatable :: z
    end function
    subroutine mmsub_( y, x)
      import :: matrix
      class( matrix), intent( out) :: y
      class( matrix), intent( in)  :: x
    end subroutine
    function mrfun_( x, y) result( z)
      import :: matrix
      class( matrix), intent( in) :: x
      real(8), intent( in)        :: y
      class( matrix), allocatable :: z
    end function
    function mcfun_( x, y) result( z)
      import :: matrix
      class( matrix), intent( in) :: x
      complex(8), intent( in)     :: y
      class( matrix), allocatable :: z
    end function
  end interface

  type, extends( matrix), public :: real_matrix
    real(8), allocatable    :: m(:,:,:)
    contains
      procedure :: plot => real_matrix_plot
      procedure :: get => real_matrix_get
      procedure :: ass => real_matrix_ass
      procedure :: add => real_matrix_add
      procedure :: sub => real_matrix_sub
      procedure :: mmmul => real_matrix_mmmul
      procedure :: mrmul => real_matrix_mrmul
      procedure :: mcmul => real_matrix_mcmul
  end type
  ! constructor
  interface real_matrix
    module procedure :: real_matrix_from_rarray
  end interface

  type, extends( matrix), public :: complex_matrix
    complex(8), allocatable :: m(:,:,:)
    contains
      procedure :: plot => complex_matrix_plot
      procedure :: get => complex_matrix_get
      procedure :: ass => complex_matrix_ass
      procedure :: add => complex_matrix_add
      procedure :: sub => complex_matrix_sub
      procedure :: mmmul => complex_matrix_mmmul
      procedure :: mrmul => complex_matrix_mrmul
      procedure :: mcmul => complex_matrix_mcmul
  end type
  ! constructor
  interface complex_matrix
    module procedure :: complex_matrix_from_rarray, complex_matrix_from_carray
  end interface

  contains

    !******************************************
    ! UTILITY FUNCTIONS
    !******************************************
    subroutine clear( M)
      class( matrix) :: M
      select type( M)
        type is( real_matrix)
          if( allocated( M%m)) deallocate( M%m)
        type is( complex_matrix)
          if( allocated( M%m)) deallocate( M%m)
      end select
      return
    end subroutine clear

    subroutine error( proc, msg)
      character(*), intent( in) :: proc, msg
      write(*,*)
      write(*,'("Error (mod_manopt_matrices/",a,"): ",a)') trim(proc), trim(msg)
      stop
    end subroutine error

    !******************************************
    ! MATRIX UTILITIES
    !******************************************
    function real_matrix_get( M) result( mat)
      class( real_matrix)  :: M
      real(8), allocatable :: mat(:,:,:)
      allocate( mat, source=M%m)
    end function

    function complex_matrix_get( M) result( mat)
      class( complex_matrix)  :: M
      complex(8), allocatable :: mat(:,:,:)
      allocate( mat, source=M%m)
    end function

    subroutine real_matrix_plot( x)
      use m_plotmat
      class( real_matrix), intent( in) :: x
      integer :: ik
      if( allocated( x%m)) then
        write(*,'("--- REAL MATRIX ---")')
        write(*,'("shape:",3i5)') shape( x%m)
        do ik = lbound( x%m, dim=3), ubound( x%m, dim=3)
          write(*,'("entry:",i5)') ik
          call plotmat( cmplx( x%m(:,:,ik), 0.d0, 8))
        end do
        write(*,'("-------------------")')
      else
        write(*,'("!!! REAL MATRIX NOT ALLOCATED !!!")')
      end if
    end subroutine

    subroutine complex_matrix_plot( x)
      use m_plotmat
      class( complex_matrix), intent( in) :: x
      integer :: ik
      if( allocated( x%m)) then
        write(*,'("--- COMPLEX MATRIX ---")')
        write(*,'("shape:",3i5)') shape( x%m)
        do ik = lbound( x%m, dim=3), ubound( x%m, dim=3)
          write(*,'("entry:",i5)') ik
          call plotmat( x%m(:,:,ik))
        end do
        write(*,'("-------------------")')
      else
        write(*,'("!!! COMPLEX MATRIX NOT ALLOCATED !!!")')
      end if
    end subroutine

    ! constructors
    function real_matrix_from_rarray( x) result( y)
      real(8), intent( in) :: x(:,:,:)
      type( real_matrix)   :: y
      if( allocated( y%m)) deallocate( y%m)
      allocate( y%m, source=x)
    end function
    function complex_matrix_from_rarray( x) result( y)
      real(8), intent( in) :: x(:,:,:)
      type( complex_matrix):: y
      if( allocated( y%m)) deallocate( y%m)
      allocate( y%m, source=cmplx( x, 0.d0, 8))
    end function
    function complex_matrix_from_carray( x) result( y)
      complex(8), intent( in) :: x(:,:,:)
      type( complex_matrix)   :: y
      if( allocated( y%m)) deallocate( y%m)
      allocate( y%m, source=x)
    end function
    ! matrix assignment
    subroutine real_matrix_ass( y, x)
      class( real_matrix), intent( out) :: y
      class( matrix), intent( in)       :: x
      select type( x)
        type is( real_matrix)
          if( allocated( y%m)) deallocate( y%m)
          allocate( y%m, source=x%m)
        class default
          call error( 'real_matrix_ass', 'Only a real matrix can be assigned to a real matrix.')
      end select
    end subroutine
    subroutine complex_matrix_ass( y, x)
      class( complex_matrix), intent( out) :: y
      class( matrix), intent( in)          :: x
      select type( x)
        type is( real_matrix)
          if( allocated( y%m)) deallocate( y%m)
          allocate( y%m, source=cmplx( x%m, 0.d0, 8))
        type is( complex_matrix)
          if( allocated( y%m)) deallocate( y%m)
          allocate( y%m, source=x%m)
        class default
          call error( 'complex_matrix_ass', 'Only a real or a complex matrix can be assigned to a complex matrix.')
      end select
    end subroutine
    ! matrix sum
    function real_matrix_add( x, y) result( z)
      class( real_matrix), intent( in) :: x
      class( matrix), intent( in)      :: y
      class( matrix), allocatable      :: z
      allocate( z, source=y)
      select type( z)
        type is( real_matrix)
          z%m = z%m + x%m
        type is( complex_matrix)
          z%m = z%m + x%m
      end select
    end function
    function complex_matrix_add( x, y) result( z)
      class( complex_matrix), intent( in) :: x
      class( matrix), intent( in)         :: y
      class( matrix), allocatable         :: z
      allocate( z, source=x)
      select type( z)
        type is( complex_matrix)
          select type( y)
            type is( real_matrix)
              z%m = z%m + y%m
            type is( complex_matrix)
              z%m = z%m + y%m
          end select
      end select
    end function
    ! matrix difference
    function real_matrix_sub( x, y) result( z)
      class( real_matrix), intent( in) :: x
      class( matrix), intent( in)      :: y
      class( matrix), allocatable      :: z
      allocate( z, source=y)
      select type( z)
        type is( real_matrix)
          z%m = x%m - z%m
        type is( complex_matrix)
          z%m = x%m - z%m
      end select
    end function
    function complex_matrix_sub( x, y) result( z)
      class( complex_matrix), intent( in) :: x
      class( matrix), intent( in)         :: y
      class( matrix), allocatable         :: z
      allocate( z, source=x)
      select type( z)
        type is( complex_matrix)
          select type( y)
            type is( real_matrix)
              z%m = z%m - y%m
            type is( complex_matrix)
              z%m = z%m - y%m
          end select
      end select
    end function
    ! matrix-matrix product (elementwise)
    function real_matrix_mmmul( x, y) result( z)
      class( real_matrix), intent( in) :: x
      class( matrix), intent( in)      :: y
      class( matrix), allocatable      :: z
      allocate( z, source=y)
      select type( z)
        type is( real_matrix)
          z%m = z%m * x%m
        type is( complex_matrix)
          z%m = z%m * x%m
      end select
    end function
    function complex_matrix_mmmul( x, y) result( z)
      class( complex_matrix), intent( in) :: x
      class( matrix), intent( in)         :: y
      class( matrix), allocatable         :: z
      allocate( z, source=x)
      select type( z)
        type is( complex_matrix)
          select type( y)
            type is( real_matrix)
              z%m = z%m * y%m
            type is( complex_matrix)
              z%m = z%m * y%m
          end select
      end select
    end function
    ! matrix-scalar product
    function real_matrix_mrmul( x, y) result( z)
      class( real_matrix), intent( in) :: x
      real(8), intent( in)             :: y
      class( matrix), allocatable      :: z
      allocate( z, source=x)
      if( abs( y - 1.d0) < 1.d-23) return
      select type( z)
        type is( real_matrix)
          z%m = y * z%m
        type is( complex_matrix)
          z%m = y * z%m
      end select
    end function
    function real_matrix_mcmul( x, y) result( z)
      class( real_matrix), intent( in) :: x
      complex(8), intent( in)          :: y
      class( matrix), allocatable      :: z
      allocate( complex_matrix :: z)
      select type( z)
        type is( complex_matrix)
          z%m = y * x%m
      end select
    end function
    function complex_matrix_mcmul( x, y) result( z)
      class( complex_matrix), intent( in) :: x
      complex(8), intent( in)             :: y
      class( matrix), allocatable         :: z
      allocate( z, source=x)
      if( abs( y - 1.d0) < 1.d-23) return
      select type( z)
        type is( complex_matrix)
          z%m = y * z%m
      end select
    end function
    function complex_matrix_mrmul( x, y) result( z)
      class( complex_matrix), intent( in) :: x
      real(8), intent( in)                :: y
      class( matrix), allocatable         :: z
      allocate( z, source=x)
      if( abs( y - 1.d0) < 1.d-23) return
      select type( z)
        type is( complex_matrix)
          z%m = y * z%m
      end select
    end function
end module mod_manopt_matrices

