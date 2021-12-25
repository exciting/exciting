subroutine diagH(n,H,S,eval,info)
      implicit none
      integer :: n,info
      complex(8) :: H(n,n),S(n,n) !,eval(n)
      real(8) :: eval(n)

      Complex (8),allocatable :: workd(:)
      real(8), allocatable :: rwork(:)
      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
      real(8) :: ta,tb

call timesec(ta)
      allocate(workd(8*n))
      allocate(rwork(6*n))
      call zhegv (1, &                ! A*x = (lambda)*B*x
                   'V', &              ! Compute eigenvalues and eigenvectors
                   'U', &              ! Upper triangle
                   n, &          ! Size of the problem
                   H, &           ! A
                   n, &          ! leading dimension of A
                   S, &           ! B
                   n, &          ! leading dimension of B
                   eval, &               ! (lambda)
                   workd, &            ! work array
                   8*n, &          ! size of lwork
                   rwork, &            ! another work array
                   info &              ! info
                   )
      deallocate(rwork,workd)
call timesec(tb)
#ifdef TIMINGS
write(*,*) 'diagH', tb-ta
#endif
!     write(*,*) info
!     write(*,*) eval
!     read(*,*)
end subroutine diagH


