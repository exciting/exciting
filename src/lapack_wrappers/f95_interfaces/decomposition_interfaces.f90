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

  end interface

end module decomposition_interfaces