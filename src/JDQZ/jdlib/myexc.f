      subroutine myexc (n, s, t, z, q, ldz, ifst, ilst)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/10/31 23:51:12 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, ldz, ifst, ilst
      double complex s(ldz,*), t(ldz,*), z(ldz,*), q(ldz,*)
c
c     .. Local ..
c
      logical tlts
      integer k, m1, m2, m3
      double complex f, f1, f2, c1, c2, r, sn
      double precision cs
c
c     .. Executable statements ..
c
      if (n.eq.1 .or. ifst.eq.ilst) return
      if (ifst.lt.ilst) then
         m1 = 0
         m2 = -1
         m3 = 1
      else
         m1 = -1
         m2 = 0
         m3 = -1
      endif
      do k = ifst+m1, ilst+m2, m3
         f = max(abs(t(k+1,k+1)),abs(s(k+1,k+1)))
         f1 = t(k+1,k+1)/f
         f2 = s(k+1,k+1)/f
         tlts = .true.
         if (abs(f1).gt.abs(f2)) tlts = .false.
         c1 = f1*s(k,k) - f2*t(k,k)
         c2 = f1*s(k,k+1) - f2*t(k,k+1)
         call zlartg (c2, c1, cs, sn, r)
         call zrot (k+1, s(1,k+1), 1, s(1,k), 1, cs, sn)
         call zrot (k+1, t(1,k+1), 1, t(1,k), 1, cs, sn)
         call zrot (n, q(1,k+1), 1, q(1,k), 1, cs, sn)
         if (tlts) then
            c1 = s(k,k)
            c2 = s(k+1,k) 
         else
            c1 = t(k,k)
            c2 = t(k+1,k) 
         endif
         call zlartg (c1, c2, cs, sn, r)
         call zrot (n-k+1, s(k,k), ldz, s(k+1,k), ldz, cs, sn)
         call zrot (n-k+1, t(k,k), ldz, t(k+1,k), ldz, cs, sn)
         call zrot (n, z(1,k), 1, z(1,k+1), 1, cs, dconjg(sn))
      enddo
c
c     --- Return
c
      end
