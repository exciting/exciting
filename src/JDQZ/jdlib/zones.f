      subroutine zones (n, x)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/07/30 21:45:46 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n
      double complex x(*)
c
c     .. Local ..
c
      integer i
c
c     .. Executable statements ..
c
      do i=1,n
         x(i) = (1.0d0,0.0d0)
      enddo
c
c     --- Return
c
      end

