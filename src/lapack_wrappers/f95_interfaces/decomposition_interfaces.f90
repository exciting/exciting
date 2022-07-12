!> FORTRAN 95 interface for LAPACK decomposition routines.
module decomposition_interfaces
  use precision, only: dp 

  interface

    ! singular value decomposition (SVD)

    ! divide and conquer algorithm for SVD
    subroutine dgesdd(jobz, m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, iwork, info)
      import :: dp 
      implicit none
      character, intent(in)    :: jobz
      integer,   intent(in)    :: m
      integer,   intent(in)    :: n
      integer,   intent(in)    :: lda
      real(dp),  intent(inout) :: A(lda, *)
      real(dp),  intent(out)   :: S(*)
      integer,   intent(in)    :: ldu
      real(dp),  intent(out)   :: U(ldu, *)
      integer,   intent(in)    :: ldvt
      real(dp),  intent(out)   :: VT(ldvt, *)
      real(dp),  intent(out)   :: work(*)
      integer,   intent(in)    :: lwork
      integer,   intent(out)   :: iwork(*)
      integer,   intent(out)   :: info
    end subroutine dgesdd

    subroutine zgesdd(jobz, m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, iwork, info)
      import :: dp 
      implicit none
      character,    intent(in)    :: jobz
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      integer,      intent(in)    :: lda
      complex(dp),  intent(inout) :: A(lda, *)
      real(dp),     intent(out)   :: S(*)
      integer,      intent(in)    :: ldu
      complex(dp),  intent(out)   :: U(ldu, *)
      integer,      intent(in)    :: ldvt
      complex(dp),  intent(out)   :: VT(ldvt, *)
      complex(dp),  intent(out)   :: work(*)
      integer,      intent(in)    :: lwork
      real(dp),     intent(out)   :: rwork(*)
      integer,      intent(out)   :: iwork(*)
      integer,      intent(out)   :: info
    end subroutine zgesdd


    ! LU factorization

    ! full pivot
    subroutine dgetc2(n, A, lda, ipiv, jpiv, info)
      import :: dp
      implicit none
      integer,  intent(in)    :: n
      integer,  intent(in)    :: lda
      real(dp), intent(inout) :: A(lda, *)
      integer,  intent(out)   :: ipiv(*)
      integer,  intent(out)   :: jpiv(*)
      integer,  intent(out)   :: info
    end subroutine dgetc2

    subroutine zgetc2(n, A, lda, ipiv, jpiv, info)
      import :: dp
      implicit none
      integer,     intent(in)    :: n
      integer,     intent(in)    :: lda
      complex(dp), intent(inout) :: A(lda, *)
      integer,     intent(out)   :: ipiv(*)
      integer,     intent(out)   :: jpiv(*)
      integer,     intent(out)   :: info
    end subroutine zgetc2

    ! row pivot
    subroutine dgetrf2(m, n, A, lda, ipiv, info)	
      import :: dp 
      implicit none
      integer, intent(in)     :: m
		  integer, intent(in)     :: n
      integer, intent(in)     :: lda
		  real(dp), intent(inout) :: A(lda, *)
		  integer, intent(out)     :: ipiv(*)
		  integer, intent(out)     :: info
	  end subroutine dgetrf2

    subroutine zgetrf2(m, n, A, lda, ipiv, info)
      import :: dp
      implicit none
      integer, intent(in)        :: m
		  integer, intent(in)        :: n
      integer, intent(in)        :: lda
		  complex(dp), intent(inout) :: A(lda, *)
		  integer, intent(out)       :: ipiv(*)
		  integer, intent(out)       :: info
	  end subroutine zgetrf2

    ! compute the inverse matrix from the LU factorization computed by *getrf
    subroutine dgetri(n, A, lda, ipiv, work, lwork, info)
      import :: dp
      implicit none
      integer,  intent(in)    :: n
      integer,  intent(in)    :: lda
      real(dp), intent(inout) :: A(lda, *)
      integer,  intent(in)    :: ipiv(*)
      real(dp), intent(out)   :: work(*)
      integer,  intent(in)    :: lwork
      integer,  intent(out)   :: info
    end subroutine dgetri

    subroutine zgetri(n, A, lda, ipiv, work, lwork, info)
      import :: dp
      implicit none
      integer,     intent(in)    :: n
      integer,     intent(in)    :: lda
      complex(dp), intent(inout) :: A(lda, *)
      integer,     intent(in)    :: ipiv(*)
      complex(dp), intent(out)   :: work(*)
      integer,     intent(in)    :: lwork
      integer,     intent(out)   :: info
    end subroutine zgetri

    ! QR factorization

    ! column pivoting
    subroutine dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
      import :: dp
      implicit none
      integer,  intent(in)    :: m
      integer,  intent(in)    :: n
      integer,  intent(in)    :: lda
      real(dp), intent(inout) :: A(lda, *)
      integer,  intent(inout) :: jpvt(*)
      real(dp), intent(out)   :: tau(*)
      real(dp), intent(out)   :: work(*)
      integer,  intent(in)    :: lwork
      integer,  intent(out)   :: info
    end subroutine

    subroutine zgeqp3(m, n, A, lda, jpvt, tau, work, lwork, rwork, info)
      import :: dp
      implicit none
      integer,     intent(in)    :: m
      integer,     intent(in)    :: n
      integer,     intent(in)    :: lda
      complex(dp), intent(inout) :: A(lda, *)
      integer,     intent(inout) :: jpvt(*)
      complex(dp), intent(out)   :: tau(*)
      complex(dp), intent(out)   :: work(*)
      integer,     intent(in)    :: lwork
      real(dp),    intent(out)   :: rwork(*)
      integer,     intent(out)   :: info
    end subroutine

    ! generate Q from elementary reflectors
    subroutine dorgqr(m, n, k, A, lda, tau, work, lwork, info)
      import :: dp
      implicit none
      integer , intent(in)    :: m
      integer , intent(in)    :: n
      integer , intent(in)    :: k
      integer , intent(in)    :: lda
      real(dp), intent(inout) :: A(lda, *)
      real(dp), intent(in)    :: tau(*)
      real(dp), intent(out)   :: work(*)
      integer , intent(in)    :: lwork
      integer , intent(out)   :: info
    end subroutine

    subroutine zungqr(m, n, k, A, lda, tau, work, lwork, info)
      import :: dp
      implicit none
      integer , intent(in)       :: m
      integer , intent(in)       :: n
      integer , intent(in)       :: k
      integer , intent(in)       :: lda
      complex(dp), intent(inout) :: A(lda, *)
      complex(dp), intent(in)    :: tau(*)
      complex(dp), intent(out)   :: work(*)
      integer , intent(in)       :: lwork
      integer , intent(out)      :: info
    end subroutine


  end interface

end module decomposition_interfaces