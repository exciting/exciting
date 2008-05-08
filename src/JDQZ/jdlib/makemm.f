      subroutine makemm (n, k, w, v, m, zm, ldm)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:33:20 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, k, ldm
      double complex w(n,*), v(n,*), m(ldm,*), zm(ldm,*)
c
c     .. Local ..
c
      integer i, j
      double complex zdotc
c
c     .. Executable statements ..
c
      do i=1,k
         do j=1,k
            if (i.eq.k.or.j.eq.k)
     $           m(i,j) = zdotc (n, w(1,i), 1, v(1,j), 1)
            zm(i,j) = m(i,j)
         enddo
      enddo
c
c     --- Return
c
      end
