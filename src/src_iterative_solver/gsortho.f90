subroutine GSortho(ld,m,n,x)
      use modmpi
      use mod_gkvector
      implicit none
      integer :: ld,n,m
      complex(8) :: x(ld,n)

      complex(8), allocatable :: zproj(:),zm(:,:),x2(:,:)
      complex(8) :: norm2,zdotc,inp,zsum,zsum2
      external :: zdotc
      real(8) :: norm
      integer :: i,j
      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
      real(8) :: ta,tb,tc,td
      complex(8), allocatable :: A(:,:), update(:,:),normlist(:),zbuf(:)
      integer, allocatable :: offset(:),ngklist(:),ibuf(:)
      integer :: lwork,lrwork,liwork,info
      complex(8), allocatable :: work(:)
      real(8), allocatable :: w(:),rwork(:)
      integer, allocatable :: iwork(:)

call timesec(ta)
if (m.gt.0) then
      Allocate(A(m,n-m))
!      Allocate(update(m,n-m))
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  m, &          ! M ... rows of op( A ) = rows of C
                  n-m, &           ! N ... cols of op( B ) = cols of C
                  ld, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  x(1,1), &           ! A
                  ld,&           ! LDA ... leading dimension of A
                  x(1,m+1), &           ! B
                  ld, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  A, &  ! C
                  m &      ! LDC ... leading dimension of C
                 )
!write(*,*) 'zgemm done'
#ifdef MPI
!      update=zero
!      call MPI_ALLREDUCE(A, update, (n-m)*m,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!      A=update
#endif
!      deallocate(update)
      allocate(update(ld,n-m))

      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  ld, &          ! M ... rows of op( A ) = rows of C
                  n-m, &           ! N ... cols of op( B ) = cols of C
                  m, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  x(1,1), & ! A
                  ld,           & ! LDA ... leading dimension of A
                  A, &           ! B
                  m,&           ! LDB ... leading dimension of B
                  zero, &          ! beta
                  update, &  ! C
                  ld &      ! LDC ... leading dimension of C
                 )
      x(1:ld,m+1:n)=x(1:ld,m+1:n)-update(1:ld,1:n-m)

!test orthogonality
#ifdef TEST_GSORTHO
      deallocate(update)
      Allocate(update(m,n-m))

      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  m, &          ! M ... rows of op( A ) = rows of C
                  n-m, &           ! N ... cols of op( B ) = cols of C
                  ld, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  x(1,1), &           ! A
                  ld,&           ! LDA ... leading dimension of A
                  x(1,m+1), &           ! B
                  ld, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  A, &  ! C
                  m &      ! LDC ... leading dimension of C
                 )
#ifdef MPI
!      update=zero
!      call MPI_ALLREDUCE(A, update, (n-m)*m,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!      A=update
#endif
      A=abs(A)
      write(*,*) 'GSortho test1',sum(A)
!      call barrier
#endif
      deallocate(A,update)

endif


if (.true.) then
      allocate(normlist(n-m))
      allocate(zbuf(n-m))

      do i=1,n-m
        normlist(i)=zdotc(ld,x(1,m+i),1,x(1,m+i),1)
      enddo
#ifdef MPI
!      zbuf=zero
!      call MPI_ALLREDUCE(normlist, zbuf, n-m, MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!      normlist=zbuf
#endif
!write(*,*) 'normlist'
      do i=1,n-m
        x(1:ld,m+i)=x(1:ld,m+i)/sqrt(dble(normlist(i)))
!write(*,*) normlist(i)
      enddo
      deallocate(normlist,zbuf)
!write(*,*) '%%%%%%%%'

      Allocate(A(n-m,n-m))
      Allocate(update(n-m,n-m))
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n-m, &          ! M ... rows of op( A ) = rows of C
                  n-m, &           ! N ... cols of op( B ) = cols of C
                  ld, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  x(1,m+1), &           ! A
                  ld,&           ! LDA ... leading dimension of A
                  x(1,m+1), &           ! B
                  ld, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  A, &  ! C
                  n-m &      ! LDC ... leading dimension of C
                 )
#ifdef MPI
!      update=zero
!      call MPI_ALLREDUCE(A, update, (n-m)*(n-m),MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!      A=update
#endif
      allocate(w(n-m))
!if (rank.eq.0) then
!write(*,*) 'A'
!do i=1,n-m
! do j=1,n-m
!   write(*,*) dble(A(j,i))
! enddo
! write(*,*)
!enddo
!endif
      allocate(work(1))
      allocate(rwork(1))
      allocate(iwork(1))
      lwork=-1
      lrwork=-1
      liwork=-1
      call zheevd ('V', 'U', n-m, A, n-m, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
      lwork=int(work(1))
      lrwork=int(rwork(1))
      liwork=iwork(1)
      deallocate(work,rwork,iwork)
      allocate(work(lwork))
      allocate(rwork(lrwork))
      allocate(iwork(liwork))
      call zheevd ('V', 'U', n-m, A, n-m, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
      deallocate(work,rwork,iwork)
!write(*,*) 'info',info
!write(*,*) 'eigval'
!write(*,*) W
!write(*,*) '------'
!write(*,*) 'eigvec'
!do i=1,n-m
! do j=1,n-m
!   write(*,*) dble(A(j,i))
! enddo
! write(*,*)
!enddo
!write(*,*) '------'

#ifdef MPI
!      call MPI_BCAST(A, (n-m)*(n-m), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr)
!      call MPI_BCAST(w, (n-m), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
#endif
      
      deallocate(update)
      allocate(update(ld,n-m))
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  ld, &          ! M ... rows of op( A ) = rows of C
                  n-m, &           ! N ... cols of op( B ) = cols of C
                  n-m, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  x(1,m+1), & ! A
                  ld,           & ! LDA ... leading dimension of A
                  A, &           ! B
                  n-m,&           ! LDB ... leading dimension of B
                  zero, &          ! beta
                  update, &  ! C
                  ld &      ! LDC ... leading dimension of C
                 )
#ifdef aggressive
    do i=1,n-m
      x(1:ld,m+i)=update(1:ld,i)/sqrt(w(i))
    enddo
#else

    x(1:ld,m+1:n)=update(1:ld,1:n-m)
      allocate(normlist(n-m))
      allocate(zbuf(n-m))

      do i=1,n-m
        normlist(i)=zdotc(ld,x(1,m+i),1,x(1,m+i),1)
      enddo
#ifdef MPI
!      zbuf=zero
!      call MPI_ALLREDUCE(normlist, zbuf, n-m, MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!      normlist=zbuf
#endif
      do i=1,n-m
        x(1:ld,m+i)=x(1:ld,m+i)/sqrt(dble(normlist(i)))
      enddo
      deallocate(normlist,zbuf)
#endif

deallocate(A,update)

endif


#ifdef TEST_GSORTHO
      Allocate(update(n,n))
      Allocate(A(n,n))
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  ld, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  x(1,1), &           ! A
                  ld,&           ! LDA ... leading dimension of A
                  x(1,1), &           ! B
                  ld, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  A, &  ! C
                  n &      ! LDC ... leading dimension of C
                 )
#ifdef MPI
!      update=zero
!      call MPI_ALLREDUCE(A, update, n*n,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!     A=update
#endif
      A=abs(A)
      write(*,*) 'GSortho test2',n,sum(A)
!      call barrier
     deallocate(A,update,w)
#endif

!zsum=sum(x(1:ld,1:n))
!call MPI_ALLREDUCE(zsum, zsum2, 1,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!write(*,*) 'GSortho done',zsum2
!call barrier



call timesec(tb)
#ifdef TIMINGS
write(*,*) 'GSortho', tb-ta
#endif
!read(*,*)
end subroutine GSortho


