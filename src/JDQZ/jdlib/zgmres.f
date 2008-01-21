      subroutine zgmres (n, x, r, mxm, eps, mxmv, 
     $     alpha, beta, k, kz, q, invqkz, ldqkz, ipiv, f, v, tp)
c
c     Programmer: Diederik R. Fokkema
c
c
      implicit none
c
c     .. Parameters ..
c
      integer mxm, mxmv, n, k, ldqkz, ipiv(*)
      double precision eps
      double complex x(*), r(*), kz(n,*), q(n,*), v(n,*),
     $     alpha, beta, invqkz(ldqkz,*), f(*), tp(*)

ctex@ \begin{manpage}{ZGMRES} 
ctex@
ctex@ \subtitle{ZGMRES} 
ctex@    ZGMRES -- Generalized Minimal Residual
ctex@    iterative method for solving linear systems $\beta A-\alpha B = -r$.
ctex@    This subroutine in specilized for use in JDQZ.
ctex@ 
ctex@ \subtitle{Declaration}
ctex@ \function{subroutine zgmres (n, x, r, mxm, eps, mxmv, a, ka, b, kb,
ctex@   alpha, beta, k, kz, mqkz, zmqkz, ldvs, q,
ctex@   lu, klu, dlu, v)}
ctex@
ctex@ \subtitle{Parameters}
ctex@    \variable{integer n} 
ctex@       On entry, n specifies the dimension of the matrix A.
ctex@       
ctex@    \variable{x} 
ctex@       double complex array of size n. 
ctex@       On exit, x is overwritten by the approximate solution.
ctex@
ctex@    \variable{r}
ctex@       double complex array of size n. Before entry, the array r 
ctex@       must contain the righthandside of the linear problem Ax=r. 
ctex@       Changed on exit.
ctex@ 
ctex@    \variable{integer mxm} 
ctex@       On entry, mxm specifies the degree of the Minimal Residual
ctex@       polynomial.
ctex@
ctex@    \variable{{double precision} eps}
ctex@       On entry, eps specifies the stop tolerance. On exit, eps contains
ctex@       the relative norm of the last residual.
ctex@
ctex@    \variable{integer mxmv}
ctex@       On Entry, mxmv specifies the maximum number of matrix 
ctex@       multiplications. On exit, mxmv contains the number of matrix
ctex@       multiplications performed by the method.
ctex@
ctex@    \variable{{double complex} zalpha}
ctex@       On entry, zalpha specifies $\alpha$. Unchanged on exit.
ctex@ 
ctex@    \variable{{double complex} zbeta}
ctex@       On entry, zbeta specifies $\beta$. Unchanged on exit.
ctex@ 
ctex@    \variable{integer k}
ctex@       On entry, k specifies the number of columns of the arrays
ctex@       kz and q.
ctex@
ctex@    \variable{z}
ctex@       double complex array z, of size (n,k). On entry the array z
ctex@       must contain the preconditioned matrix Z.
ctex@
ctex@    \variable{mqkz}
ctex@       double complex array mqkz, of size (ldvs,k). On entry the array 
ctex@       mqkz must contain the matrix Q'*KZ.
ctex@
ctex@    \variable{zmqkz}
ctex@       double complex array zmqkz, of size (ldvs,k). Workspace. Used to
ctex@       copy mqkz.
ctex@
ctex@    \variable{q}
ctex@       double complex array q, of size (n,k). On entry the array q
ctex@       must contain the preconditioned matrix Q.
ctex@
ctex@    \variable{v}
ctex@       double complex array of size (n,mxm+1). Workspace.
ctex@
ctex@ \subtitle{Description}
ctex@    ***
ctex@
ctex@ \subtitle{See Also}
ctex@    ***
ctex@
ctex@ \subtitle{References}
ctex@    ***
ctex@ 
ctex@ \subtitle{Bugs}
ctex@    ***
ctex@
ctex@ \subtitle{Author}
ctex@     Diederik R.\ Fokkema
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
ctex@
c
c     .. Local ..
c
      logical restrt, loop
      integer maxm
      parameter (maxm = 100)

      integer i, m, m1, nmv
      double precision rnrm0, rnrm, eps1, c(maxm)
      double complex hh(maxm,maxm-1), rs(maxm), s(maxm), y(maxm), rcs
      double complex zero, one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      double precision  dznrm2
      double complex ztmp, zdotc
c
c     .. Executable Statements ..
c
      if ( mxm .gt. maxm-1 )
     $     call error ('mxm larger than maxm (zgmres)')
c
c     --- Initialize first residue
c
      call zzeros (n, x)
      call zscal (n, -one, r, 1)
      call psolve (n, r, k, q, kz, invqkz, ldqkz,
     $     ipiv, f)
c
c     --- initialize loop
c
      rnrm0  = dznrm2 (n,r,1)
      rnrm = rnrm0
      eps1  = eps * rnrm

      nmv = 0

      call zcopy (n, r, 1, v(1,1), 1)
c         
c     --- outer loop
c
 1000 restrt = ( nmv .lt. mxmv .and. rnrm .gt. eps1 )
      if ( restrt ) then
         ztmp = one / rnrm
         call zscal (n, ztmp, v(1,1), 1)
         rs(1) = rnrm
c
c     --- inner loop
c
         m = 0
 2000    loop = (nmv.lt.mxmv .and. m.lt.mxm .and. rnrm.gt.eps1)
         if (loop) then
            m  = m + 1
            m1 = m + 1
            call jdqzmv (n, v(1,m), v(1,m1), tp, alpha, beta)
            call psolve (n, v(1,m1), k, q, kz, invqkz,
     $           ldqkz, ipiv, f)
            nmv = nmv+1 
            do i = 1,m
               ztmp = zdotc (n, v(1,i), 1, v(1,m1), 1)
               hh(i,m) = ztmp
               call zaxpy (n, (-ztmp), v(1,i), 1, v(1,m1), 1)
            enddo
            ztmp = dznrm2( n, v(1,m1), 1 )
            hh(m1,m) = ztmp
            call zscal ( n, (one/ztmp), v(1,m1), 1 )
            do i = 1,m-1
               call zrot (1, hh(i,m), 1, hh(i+1,m), 1, c(i), s(i))
            enddo
            call zlartg (hh(m,m), hh(m1,m), c(m), s(m), rcs )
            hh(m,m) = rcs
            hh(m1,m) = zero
            rs(m1) = zero
            call zrot (1, rs(m), 1, rs(m1), 1, c(m), s(m))
            rnrm = abs (rs(m1))
            goto 2000
         endif
c
c     --- compute approximate solution x
c
         call zcopy ( m, rs, 1, y, 1 )
         call ztrsv ( 'u', 'n', 'n', m, hh, maxm, y, 1 )
         call zgemv ( 'n', n, m, one, v, n, y, 1, one, x, 1 )
c
c     --- compute residual for restart
c
         call jdqzmv (n, x, v(1,2), tp, alpha, beta)
         call psolve (n, v(1,2), k, q, kz, invqkz,
     $        ldqkz, ipiv, f)
         call zcopy (n, r, 1, v(1,1), 1)
         call zaxpy (n, -one, v(1,2), 1, v(1,1), 1)
         rnrm = dznrm2 (n, v(1,1), 1)

         goto 1000
      endif
c
c     --- return
c
      eps = rnrm/rnrm0
      mxmv = nmv

      return
      end
