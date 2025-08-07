module linear_system_interfaces
  use precision, only: dp, i32

  implicit none

  private

  public :: zgelsd, zposv

  interface

    subroutine zgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, &
                      rwork, iwork, info)
      import :: dp, i32
      integer(i32), intent(in)    :: m
      integer(i32), intent(in)    :: n
      integer(i32), intent(in)    :: nrhs
      complex(dp), intent(inout)  :: a(lda, *)
      integer(i32), intent(in)    :: lda
      complex(dp), intent(inout)  :: b(ldb, *)
      integer(i32), intent(in)    :: ldb
      real(dp), intent(out)       :: s(*)
      real(dp), intent(in)        :: rcond
      integer(i32), intent(out)   :: rank
      complex(dp), intent(out)    :: work(*)
      integer(i32), intent(in)    :: lwork
      real(dp), intent(out)       :: rwork(*)
      integer(i32), intent(out)   :: iwork(*)
      integer(i32), intent(out)   :: info
    end subroutine zgelsd
    subroutine zposv( uplo, n, nrhs, A, lda, B, ldb, info )
      import :: dp, i32
      character, intent(in)       :: uplo
      integer(i32), intent(in)    :: n
      integer(i32), intent(in)    :: nrhs
      complex(dp), intent(inout)  :: A( lda, * )
      integer(i32), intent(in)    :: lda
      complex(dp), intent(inout)  :: B( lda, * )
      integer(i32), intent(in)    :: ldb
      integer(i32), intent(out)   :: info
    end subroutine
  end interface
end module