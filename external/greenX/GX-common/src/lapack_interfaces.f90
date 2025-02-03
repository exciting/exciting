! **************************************************************************************************
!  Copyright (C) 2020-2023 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief Declare interfaces for BLAS/LAPACK subroutines
! **************************************************************************************************
module lapack_interfaces
  implicit none


  interface blas3
     subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
       use kinds, only: dp
       implicit none
       integer                      :: len
       character(len=1), intent(in) :: transa, transb
       integer, intent(in)          :: m, n, k, lda, ldb, ldc
       real(kind=dp), intent(in)    :: alpha, beta
       real(kind=dp), intent(in)    :: a(lda, *), b(ldb, *)
       real(kind=dp), intent(inout) :: c(ldc, *)
     end subroutine dgemm
  end interface

  interface svd
     subroutine dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
       use kinds, only: dp
       implicit none
       integer                      :: len
       character(len=1), intent(in) :: jobz
       integer, intent(in)          :: m, n, lda, ldu, ldvt, lwork
       integer, intent(out)         :: info
       real(kind=dp), intent(out)   :: s(*), u(ldu, *), vt(ldvt, *), work(*)
       real(kind=dp), intent(inout) :: a(lda, *)
       integer, intent(inout)       :: iwork(*)
     end subroutine dgesdd
  end interface

  !interface num_prec
  !   real(kind=dp) function dlamch(cmach)
  !     implicit none             
  !     character(len=1), intent(in) :: cmach
  !   end function dlamch 
  !end interface

  interface diag
     subroutine dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &
                            m, w, z, ldz, work, lwork, iwork, ifail, info) 
       use kinds, only: dp     
       implicit none
       character(len=1), intent(in) :: jobz, range, uplo
       integer, intent(in)          :: il, iu, lda, ldz, lwork, m, n
       integer, intent(out)         :: info
       real(kind=dp), intent(in)    :: abstol, vl, vu

       integer, intent(in)          :: ifail(*), iwork(*)
       real(kind=dp), intent(inout) :: a(lda, *), w(*), work(*), z(ldz,*)
     end subroutine dsyevx
  end interface

end module lapack_interfaces 
