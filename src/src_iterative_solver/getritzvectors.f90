subroutine getritzvectors(l,m,n,A,B,C,nstart,npw)
      use mod_gkvector
      use modmpi
      implicit none
      integer :: l,m,n,nstart,npw
      complex(8) :: A(l,m),B(m,m),C(l,n)

      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
      real(8) :: ta,tb
      integer :: nlo,i,j
      integer, allocatable :: offset(:),ngklist(:),ibuf(:)

call timesec(ta)
if (.true.) then
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  l, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  m, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  A, &           ! A
                  l,&           ! LDA ... leading dimension of A
                  B(1,nstart), &           ! B
                  m, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  C, &  ! C
                  l &      ! LDC ... leading dimension of C
                  )
else
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  npw, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  m, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  A(1,1), &           ! A
                  l,&           ! LDA ... leading dimension of A
                  B(1,nstart), &           ! B
                  m, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  C(1,1), &  ! C
                  l &      ! LDC ... leading dimension of C
                  )
      nlo=l-npw
      if (nlo.ne.0) then
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    nlo, &          ! M ... rows of op( A ) = rows of C
                    n, &           ! N ... cols of op( B ) = cols of C
                    m, &          ! K ... cols of op( A ) = rows of op( B )
                    one, &          ! alpha
                    A(npw+1,1), &           ! A
                    l,&           ! LDA ... leading dimension of A
                    B(1,nstart), &           ! B
                    m, &          ! LDB ... leading dimension of B
                    zero, &          ! beta
                    C(npw+1,1), &  ! C
                    l &      ! LDC ... leading dimension of C
                    )
     endif

endif


!endif

call timesec(tb)
#ifdef TIMINGS
write(*,*) 'getritzvectors', tb-ta
#endif
end subroutine getritzvectors


