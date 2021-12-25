subroutine innerproduct(l,m,n,A,B,C)
      use modmpi
      use mod_gkvector
      implicit none
      integer :: l,m,n
      complex(8) :: A(l,m),B(l,n),C(m,n)

      Complex(8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
      real(8) :: ta,tb
      complex(8), allocatable :: buf(:,:)


      call timesec(ta)
      C=zero
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  m, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  l, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  A(1,1), &           ! A
                  l,&           ! LDA ... leading dimension of A
                  B(1,1), &           ! B
                  l, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  C, &  ! C
                  m &      ! LDC ... leading dimension of C
                  )
#ifdef MPI
!        allocate(buf(m,n))
!        buf=zero
!        call MPI_ALLREDUCE(C, buf, m*n,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!        C=buf
!        deallocate(buf)
#endif       
!write(*,*) 'sumC2',sum(C)

     call timesec(tb)
#ifdef TIMINGS
     write(*,*) 'innerproduct', tb-ta
#endif

end subroutine innerproduct


