subroutine diagH2(n,H,S,eval,info)
      use modmpi
      implicit none
      integer :: n,info
      complex(8) :: H(n,n),S(n,n)
      real(8) :: eval(n)

      Complex (8),allocatable :: work(:)
      real(8), allocatable :: rwork(:)
      integer, allocatable :: iwork(:)
      integer :: lwork,lrwork,liwork
      complex(8) :: worksize
      real(8) :: rworksize
      integer :: iworksize,itype
      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
      real(8) :: ta,tb

      itype=1
      lwork=-1
      lrwork=-1
      liwork=-1
      call zhegvd(itype, &              ! A*x = (lambda)*B*x
                  'V', &              ! Compute eigenvalues and eigenvectors
                  'U', &              ! Upper triangle
                   n,  &              ! Size of the problem
                   H,  &
                   n,  & 
                   S,  &
                   n,  & 
                   eval, &
                   worksize,   &
                   lwork,  &
                   rworksize,  & 
                   lrwork, &
                   iworksize,  &
                   liwork, &
                   info )
      lwork=int(dble(worksize))
      lrwork=int(rworksize)
      liwork=iworksize
      allocate(work(lwork))
      allocate(rwork(lrwork))
      allocate(iwork(liwork))

      call zhegvd(itype, &              ! A*x = (lambda)*B*x
                  'V', &              ! Compute eigenvalues and eigenvectors
                  'U', &              ! Upper triangle
                   n,  &              ! Size of the problem
                   H,  &
                   n,  & 
                   S,  &
                   n,  & 
                   eval, &
                   work,   &
                   lwork,  &
                   rwork,  & 
                   lrwork, &
                   iwork,  &
                   liwork, &
                   info )
      deallocate(rwork,work,iwork)

!      write(*,*) 'info',info
end subroutine diagH2

