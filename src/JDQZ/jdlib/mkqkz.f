      subroutine mkqkz (n, k, q, kq, qkz, invqkz, ldqkz, ipiv)
c
c     Coded by Diederik Fokkema
c
c     $Id: mkqkz.f,v 1.1 1995/08/05 09:08:21 caveman Exp $
c
c     Time-stamp: <95/08/03 23:52:52 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, k, ldqkz, ipiv(*)
      double complex q(n,*), kq(n,*), qkz(ldqkz,*), invqkz(ldqkz,*)
c
c     .. local ..
c
      integer i, j, info
      double complex zdotc
c
c     .. Executable statements ..
c
      do i=1,k
         do j=1,k
            if (i.eq.k.or.j.eq.k) qkz(i,j) =
     $           zdotc (n, q(1,i), 1, kq(1,j), 1)
            invqkz(i,j) = qkz(i,j)
         enddo
      enddo
      call zgetrf (k, k, invqkz, ldqkz, ipiv, info)
c
c     --- Return
c
      end

