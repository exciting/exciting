      program zjdqz
c
      implicit none

c     ... work space

      integer lwork, n
      parameter (n = 50)
      double complex zwork(200 000), alpha(20), beta(20), 
     $               eivec(n,20), tmp(n), residu(n)

      integer kmax, jmax, jmin, maxstep, method, m, l, maxnmv, 
     $     order, testspace, j
      double precision tol, lock, dznrm2
      double complex target

c     ... Start of program

c....    Number of double complex numbers in workspace:     
      lwork = 200 000
      target = (0.5d0,0.0)
      tol = 1.d-9
      kmax = 5
      jmax = 15
      jmin = 10
      maxstep = 200
      lock = 1.d-8
c...     order =  0: nearest to target
c...     order = -1: smallest real part
c...     order =  1: largest real part
c...     order = -2: smallest complex part
c...     order =  2: largest complex part
      order = 0
c...     method = 1: gmres(m)
c...     method = 2: cgstab(l)
      method = 2
c...     for gmres(m):
      m = 10
c...     for cgstab(l):
      l= 4
c...     maximum number of matvecs in cgstab or gmres
      maxnmv = 100
      testspace = 3
c...     Testspace 1: w = "Standard Petrov" * v (Section 3.1.1)
c...     Testspace 2: w = "Standard 'variable' Petrov" * v (Section 3.1.2)
c...     Testspace 3: w = "Harmonic Petrov" * v (Section 3.5.1)
      print '(/a/,80("="))', 'input:'
      print '(a20,1p," (",e10.3,",",e10.3," )")', 'target:', target
      print '(a20,1p,e12.3)', 'tol:', tol
      print '(a20,i12)', 'kmax:', kmax
      print '(a20,i12)', 'jmax:', jmax
      print '(a20,i12)', 'jmin:', jmin
      print '(a20,i12)', 'maxstep:', maxstep
      print '(a20,1p,e12.3)', 'lock:', lock
      print '(a20,i12)', 'order:', order
      print '(a20,i12)', 'method:', method
      print '(a20,i12)', 'gmres m:', m
      print '(a20,i12)', 'bicgstab l:', l
      print '(a20,i12)', 'mxmv:', maxnmv
      print '(a20,i12)', 'testspace:', testspace

      call jdqz(alpha, beta, eivec, .true., n, target, tol, 
     $     kmax, jmax, jmin,
     $     method, m, l, maxnmv, maxstep,
     $     lock, order, testspace, zwork, lwork )

      do j = 1, kmax
         call amul  ( n, eivec(1,j), residu )
         call zscal ( n, beta(j), residu, 1)
         call bmul  ( n, eivec(1,j), tmp )
         call zaxpy( n, -alpha(j), tmp, 1, residu, 1 )
         print '("lambda(",i2,"): (",1p,e11.4,",",e11.4,
     $        " )")', j,alpha(j)/beta(j)
         print '(a30,d13.6)', '||beta Ax - alpha Bx||:',
     $          dznrm2( n, residu, 1 )
      end do
      end
c
      subroutine precon( neq, q )
      double complex q(*)
      integer neq
      end
      subroutine amul( neq, q, r )
      double complex q(*), r(*)
      integer neq, i
      do 10, i = 3, neq-1
         r(i) = q(i-2)+2d0*q(i)+q(i+1)
   10 continue
      r(1)   = 2d0*q(1)  +q(2)
      r(2)   = 2d0*q(2)  +q(3)
      r(neq) = q(neq-2)+2d0*q(neq)
      end
      subroutine bmul( neq, q, r )
      double complex q(*), r(*)
      integer neq, i
      do 10, i = 1, neq
	 r(i) = q(i)
   10 continue
      end

