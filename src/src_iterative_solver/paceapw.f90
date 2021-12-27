subroutine paceapw(n,npw,nwf,nproj,prefactor,P,x,Hx)
      use modmain, only : nmatmax,zzero,zone
      use mod_hybrids, only : vnlmat
      implicit none
      integer ::n,nwf,npw,nproj
      real(8) :: prefactor
      complex(8) :: x(n,nwf),Sx(n,nwf),Hx(n,nwf), P(nproj,nmatmax)

      complex(8), allocatable :: zax(:,:)
      complex(8) :: zp
      real(8) :: ta,tb

call timesec(ta)

if (.true.) then
      allocate(zax(nproj,nwf))
      zp=-prefactor
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  nproj, &          ! M ... rows of op( A ) = rows of C
                  nwf, &           ! N ... cols of op( B ) = cols of C
                  npw, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  P, &           ! A
                  nproj,&           ! LDA ... leading dimension of A
                  x, &           ! B
                  n, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  zax, &  ! C
                  nproj &      ! LDC ... leading dimension of C
                 )
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n, &          ! M ... rows of op( A ) = rows of C
                  nwf, &           ! N ... cols of op( B ) = cols of C
                  nproj, &          ! K ... cols of op( A ) = rows of op( B )
                  zp, &          ! alpha
                  P, &           ! A
                  nproj,&           ! LDA ... leading dimension of A
                  zax, &           ! B
                  nproj, &          ! LDB ... leading dimension of B
                  zone, &          ! beta
                  Hx, &  ! C
                  n &      ! LDC ... leading dimension of C
                )
      deallocate(zax)
else
      zp=prefactor
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n, &          ! M ... rows of op( A ) = rows of C
                  nwf, &           ! N ... cols of op( B ) = cols of C
                  n, &          ! K ... cols of op( A ) = rows of op( B )
                  zp, &          ! alpha
                  vnlmat(1,1,1), &           ! A
                  nmatmax,&           ! LDA ... leading dimension of A
                  x, &           ! B
                  n, &          ! LDB ... leading dimension of B
                  zone, &          ! beta
                  Hx, &  ! C
                  n &      ! LDC ... leading dimension of C
                  )
endif



call timesec(tb)

#ifdef TIMINGS
write(*,*) 'paceapw',tb-ta,td-tc
#endif
end subroutine paceapw

