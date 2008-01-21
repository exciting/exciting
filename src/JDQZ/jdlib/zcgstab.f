      subroutine zcgstabl (n, x, r, l, eps, mxmv,
     $     zalpha, zbeta, nk, kz, qq, invqkz, ldqkz, jpiv, f,
     $     work, lwork)
c
c     Programmer: Diederik R. Fokkema
c
c
      implicit none
c
c     .. Parameters ..
c
      integer l, n, mxmv, nk, lwork, ldqkz, jpiv(*)
      double precision  eps
      double complex x(*), r(*), kz(n,*), qq(n,*), work(n,*),
     $     zalpha, zbeta, invqkz(ldqkz,*), f(*)

c
c     .. Local ..
c
      integer mxl
      parameter (mxl = 32)

      integer i, j, k, m, ipiv(mxl), nmv, info

      integer u, q, w, rr, xp, bp
      logical rcomp, xadd

      double precision maxnorm, delta, bpnorm, tr0norm, trlnorm
      double complex varrho, hatgamma

      double precision rnrm, rnrm0, eps1, dznrm2
      double complex alpha, beta, omega, gamma, rho0, rho1,
     $     yr(mxl), ys(mxl), z(mxl,mxl), ztmp

      double complex zero, one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      double complex zdotc
c
c     .. Executable statements ..
c
      if (mxmv.eq.0) return

      if (l.gt.mxl)
     $     call error ('l exceeds mxl (zcgstabl)')

      u  = 1
      q  = u + (l+1)
      rr = q + (l+1)
      xp = rr + 1
      bp = xp + 1
      w = bp + 1

      if (w.gt.lwork)
     $     call error ('workspace too small (zcgstabl)')

c
c     --- set x to zero and compute first residue
c
      call zzeros (n, x)
      call zscal (n, -one, r, 1)
      call zcopy (n, r, 1, work(1,rr), 1)
      call psolve (n, r, nk, qq, kz, invqkz, ldqkz,
     $     jpiv, f)

c
c     --- initialize loop
c
      nmv = 0

      rnrm0 = dznrm2 (n, r, 1)
      rnrm = rnrm0
      eps1 = eps * rnrm0

      call zcopy (n, x, 1, work(1,xp), 1)
      call zcopy (n, r, 1, work(1,bp), 1)
      maxnorm = 0d0
      bpnorm = rnrm
      rcomp = .false.
      xadd = .false.
      delta = 1d-2

      m = 0
      rho0 = one
      alpha = one
      beta = zero
      omega = one

      call zzeros (n*(l+1), work(1,u))
      call zzeros (n*(l+1), work(1,q))
      call zcopy (n, r, 1, work(1,q), 1)
c
c     --- loop
c
 1000 continue
      m = m + l
c
c     --- BiCG part
c
      rho0 = -omega * rho0
      do k=1,l
         rho1 = zdotc (n, work(1,rr), 1, work(1,q+k-1), 1)
         beta = alpha * (rho1 / rho0)
         rho0 = rho1
         beta = beta
         do j=0,k-1
            call zxpay (n, work(1,q+j), 1, (-beta), work(1,u+j), 1)
         enddo
         call jdqzmv (n, work(1,u+k-1), work(1,u+k), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,u+k), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         gamma = zdotc (n, work(1,rr), 1, work(1,u+k), 1)
         alpha = rho0 / gamma
         do j=0,k-1
            call zaxpy (n, (-alpha), work(1,u+j+1), 1, work(1,q+j), 1)
         enddo
         call jdqzmv (n, work(1,q+k-1), work(1,q+k), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,q+k), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         call zaxpy (n, alpha, work(1,u), 1, x, 1)
         rnrm = dznrm2 (n, work(1,q), 1)
         maxnorm = max (maxnorm, rnrm)
         nmv = nmv+2
      enddo
c
c     --- MR part + Maintaining the convergence
c
      do i=1,l-1
         do j=1,i
            ztmp = zdotc (n, work(1,q+i), 1, work(1,q+j), 1)
            z(i,j) = ztmp
            z(j,i) = dconjg(ztmp)
         enddo
         yr(i) = zdotc (n, work(1,q+i), 1, work(1,q), 1)
         ys(i) = zdotc (n, work(1,q+i), 1, work(1,q+l), 1)
      enddo
      call zgetrf (l-1, l-1, z, mxl, ipiv, info)
      call zgetrs ('n', l-1, 1, z, mxl, ipiv, yr, mxl, info)
      call zgetrs ('n', l-1, 1, z, mxl, ipiv, ys, mxl, info)
      call zcopy (n, work(1,q), 1, r, 1)
      call zgemv ('n', n, l-1, (-one), work(1,q+1), n, yr, 1, one,
     $     r, 1)
      call zgemv ('n', n, l-1, (-one), work(1,q+1), n, ys, 1, one,
     $     work(1,q+l), 1)
      tr0norm = dznrm2 (n, r, 1)
      trlnorm = dznrm2 (n, work(1,q+l), 1)
      varrho = zdotc (n, work(1,q+l), 1, r, 1) / (tr0norm*trlnorm)
      hatgamma = varrho/abs(varrho) * max(abs(varrho),7d-1)
      hatgamma = (tr0norm/trlnorm)*hatgamma
      yr(l) = zero
      ys(l) = -one
      call zaxpy (l, (-hatgamma), ys, 1, yr, 1)

      omega = yr(l)
      call zgemv ('n', n, l, one, work(1,q), n, yr, 1, one, x, 1)
      call zgemv ('n', n, l, (-one), work(1,u+1), n, yr, 1, one,
     $     work(1,u), 1)
      call zaxpy (n, (-hatgamma), work(1,q+l), 1, r, 1)
      call zcopy (n, r, 1, work(1,q), 1)
c
c     --- reliable update
c
      rnrm = dznrm2 (n, work(1,q), 1)
      maxnorm = max (maxnorm, rnrm)
      xadd = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.maxnorm)
      rcomp = ((rnrm.lt.delta*maxnorm.and.rnrm0.lt.maxnorm).or.xadd)

      if (rcomp) then
         call jdqzmv (n, x, work(1,q), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,q), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         call zxpay (n, work(1,bp), 1, -one, work(1,q), 1)
         maxnorm = rnrm
         if (xadd) then
            call zaxpy (n, one, x, 1, work(1,xp), 1)
            call zzeros (n, x)
            call zcopy (n, work(1,q), 1, work(1,bp), 1)
            bpnorm = rnrm
         endif
      endif
         
      if (nmv.lt.mxmv .and. rnrm.gt.eps1) goto 1000

      call zaxpy (n, one, work(1,xp), 1, x, 1)
c
c     --- return
c
      mxmv = nmv
      eps = rnrm/rnrm0

      return
      end
