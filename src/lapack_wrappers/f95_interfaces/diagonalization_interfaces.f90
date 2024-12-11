!> FORTRAN 95 interface for LAPACK diagonalization routines.
module diagonalization_interfaces
  use precision, only: dp, i32

  implicit none

  private

  public :: dsyev, &
            dstedc, &
            dsygvx, &
            zheev, &
            zhegvx

  interface
      ! Symmetric matrix
      subroutine dsyev( jobz, uplo, n, A, lda, w, work, lwork, info )
        import :: dp, i32
        implicit none
        character, intent(in)     :: jobz
        character, intent(in)     :: uplo
        integer(i32), intent(in)  :: n
        real(dp), intent(inout)   :: A(lda, *)
        integer(i32), intent(in)  :: lda
        real(dp), intent(out)     :: w(*)
        real(dp), intent(out)     :: work(*)
        integer(i32), intent(in)  :: lwork
        integer(i32), intent(out) :: info
      end subroutine

      ! Tridiagonal matrix
      subroutine dstedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
        import :: dp 
        implicit none
        character, intent(in)    :: compz
        integer,   intent(in)    :: n
        real(dp),  intent(inout) :: d(*)
        real(dp),  intent(inout) :: e(*)
        integer,   intent(in)    :: ldz
        real(dp),  intent(inout) :: z(ldz, *)
        real(dp),  intent(out)   :: work(*)
        integer,   intent(in)    :: lwork
        integer,   intent(out)   :: iwork(*)
        integer,   intent(in)    :: liwork
        integer,   intent(out)   :: info
      end subroutine dstedc

      ! Real generalized symmetric-definite eigenvalue problem
      subroutine dsygvx(itype, jobz, range, uplo, n, A, lda, B, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
        import :: dp
        implicit none
        integer,   intent(in)    :: itype
        character, intent(in)    :: jobz
        character, intent(in)    :: range
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n
        real(dp),  intent(inout) :: A(lda, *)
        integer,   intent(in)    :: lda
        real(dp),  intent(inout) :: B(ldb, *)
        integer,   intent(in)    :: ldb
        real(dp),  intent(in)    :: vl
        real(dp),  intent(in)    :: vu
        integer,   intent(in)    :: il
        integer,   intent(in)    :: iu
        real(dp),  intent(in)    :: abstol
        integer,   intent(out)   :: m
        real(dp),  intent(out)   :: w(*)
        real(dp),  intent(out)   :: z(ldz, *)
        integer,   intent(in)    :: ldz
        real(dp),  intent(out)   :: work(*)
        integer,   intent(in)    :: lwork
        integer,   intent(out)   :: iwork(*)
        integer,   intent(out)   :: ifail(*)
        integer,   intent(out)   :: info
      end subroutine dsygvx

      ! Hermitian matrix
      subroutine zheev( jobz, uplo, n, A, lda, w, work, lwork, rwork, info )
        import :: dp, i32
        implicit none
        character, intent(in)      :: jobz
        character, intent(in)      :: uplo
        integer(i32), intent(in)   :: n
        complex(dp), intent(inout) :: A(lda, *)
        integer(i32), intent(in)   :: lda
        real(dp), intent(out)      :: w(*)
        complex(dp), intent(out)   :: work(*)
        integer(i32), intent(in)   :: lwork
        real(dp), intent(out)      :: rwork(*)
        integer(i32), intent(out)  :: info
      end subroutine

      ! Complex generalized Hermitian-definite eigenvalue problem
      subroutine zhegvx(itype, jobz, range, uplo, n, A, lda, B, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
        import :: dp, i32
        implicit none
        integer(i32), intent(in)   :: itype
        character, intent(in)      :: jobz
        character, intent(in)      :: range
        character, intent(in)      :: uplo
        integer(i32), intent(in)   :: n
        complex(dp), intent(inout) :: A(lda, *)
        integer(i32), intent(in)   :: lda
        complex(dp), intent(inout) :: B(ldb, *)
        integer(i32), intent(in)   :: ldb
        real(dp), intent(in)       :: vl
        real(dp), intent(in)       :: vu
        integer(i32), intent(in)   :: il
        integer(i32), intent(in)   :: iu
        real(dp), intent(in)       :: abstol
        integer(i32), intent(out)  :: m
        real(dp), intent(out)      :: w(*)
        complex(dp), intent(out)   :: z(ldz, *)
        integer(i32), intent(in)   :: ldz
        complex(dp), intent(out)   :: work(*)
        integer(i32), intent(in)   :: lwork
        real(dp), intent(out)      :: rwork(*)
        integer(i32), intent(out)  :: iwork(*)
        integer(i32), intent(out)  :: ifail(*)
        integer(i32), intent(out)  :: info
      end subroutine zhegvx

  end interface

end module diagonalization_interfaces