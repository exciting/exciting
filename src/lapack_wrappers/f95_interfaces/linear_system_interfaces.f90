module linear_system_interfaces
  use precision, only: dp, i32

  implicit none

  private

  public :: zposv

  interface
    subroutine zposv( uplo, n, nrhs, A, lda, B, ldb, info )
      import :: dp, i32
      character, intent(in)       :: uplo
      integer(i32), intent(in)    :: n
      integer(i32), intent(in)    :: nrhs
      complex(dp), intent(inout)  :: A( lda, * )
      integer(i32), intent(in)    ::	lda
      complex(dp), intent(inout)  :: B( lda, * )
      integer(i32), intent(in)    :: ldb
      integer(i32), intent(out)   :: info
    end subroutine
  end interface
end module