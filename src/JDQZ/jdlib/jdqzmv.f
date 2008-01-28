      subroutine jdqzmv (n, x, y, work, alpha, beta)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n
      double complex alpha, beta, x(*), y(*), work(*) 
c
c     .. Local ..
c
      integer i
c     .. Executable statements ..
c
      call amul( n, x, work )
      call bmul( n, x, y )
      do i = 1, n
         y(i) = beta*work(i) - alpha*y(i)
      end do
c
      end

