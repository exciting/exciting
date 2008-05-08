      subroutine qzsort (ta, tb, k, s, t, z, q, ldz, alpha, beta,
     $     order)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:34:03 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer k, ldz, order
      double complex ta, tb, s(ldz,*), t(ldz,*), z(ldz,*),
     $     q(ldz,*), alpha(*), beta(*)
c
c     .. Local ..
c
      integer i, j, select
c
c     .. Executable statements ..
c
      do i = 1,k
         do j = 1,k
            alpha(j) = s(j,j)
            beta(j) = t(j,j)
         enddo
         j = select (k-i+1, ta, tb, alpha(i), beta(i), order) + i-1
         call myexc (k, s, t, z, q, ldz, j, i)
      enddo
c
c     --- Return
c
      end
