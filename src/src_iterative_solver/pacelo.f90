subroutine pacelo(n,npw,nwf,nproj,prefactor,P,x,Hx)
      use modmain, only : nmatmax,zzero,zone
      implicit none
      integer ::n,nwf,npw,nproj
      real(8) :: prefactor
      complex(8) :: x(n,nwf),Sx(n,nwf),Hx(n,nwf), P(nproj,nmatmax)

      integer :: nlo
      complex(8), allocatable :: zax(:,:)
      complex(8) :: zp
      real(8) :: ta,tb

call timesec(ta)
      if (n.ne.npw) then
        nlo=n-npw
        zp=-prefactor
        allocate(zax(nproj,nwf))
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    nproj, &          ! M ... rows of op( A ) = rows of C
                    nwf, &           ! N ... cols of op( B ) = cols of C
                    nlo, &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    P(1,npw+1), &           ! A
                    nproj,&           ! LDA ... leading dimension of A
                    x(npw+1,1), &           ! B
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
      endif
call timesec(tb)

#ifdef TIMINGS
write(*,*) 'paceapw',tb-ta,td-tc
#endif
end subroutine pacelo

