      subroutine zmgs (n, k, v, w, job )
c
c     Coded by Diederik Fokkema
c     Modified 05-23-96: M. Kooper: job =1 corrected, array YWORK added
c
c     .. Parameters ..
c
      implicit none
      integer n, k, job
      double complex v(n,*), w(*)
c
c     .. Local ..
c
      integer i
      double precision s0, s1, dznrm2
      double complex znrm
c
c     .. Executable statements ..
c
      s1 = dznrm2 (n, w, 1)
      do i=1, k
	 s0 = s1
	 call zortho (n, v(1,i), w, s0, s1, znrm)
      enddo
      if (job.eq.0) then
	 return
      else
         znrm  = 1.d0/s1
         call zscal (n, znrm, w, 1)
      endif
c
c     --- Return
c
      return
      end

      subroutine zortho (n, v, w, s0, s1, znrm)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n
      double precision s0, s1
      double complex v(*), w(*), znrm
c
c     .. Local ..
c
      double precision kappa, dznrm2
      double complex ztmp, zdotc
      parameter (kappa=1d2)
c
c     .. Executable statements ..
c
      znrm = zdotc (n, v, 1, w, 1)
      call zaxpy (n, (-znrm), v, 1, w, 1)
      s1 = dznrm2 (n, w, 1)
      if (s1.gt.s0/kappa) then
	 goto 100
      else
	 s0 = s1
	 ztmp = zdotc (n, v, 1, w, 1)
         znrm = znrm + ztmp
	 call zaxpy (n, (-ztmp), v, 1, w, 1)
	 s1 = dznrm2 (n, w, 1)
	 if (s1.gt.s0/kappa) then
	    goto 100
	 else
	    call error ('zero vector in zmgs')
	 endif
      endif
 100  continue
c
c     --- Return
c
      return
      end

