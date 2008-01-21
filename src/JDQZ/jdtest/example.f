      program example
c
      integer          lwork, n
      parameter        (n = 100)
      double complex   zwork(200 000), alpha(20), beta(20), 
     $                 eivec(n,20), tmp(n), residu(n)
c
      integer          kmax, jmax, jmin, maxstep, method, m, l, maxnmv, 
     $                 order, testspace, j
      double precision tol, lock, dznrm2
      logical          wanted
      double complex   target
      real             elapse
      real             etime, tarray(2)
c
      target = (31.,0.0)
      tol = 1.d-9
      kmax = 5
      jmax = 20
      jmin = 10
      maxstep = 1000
      lock = 1.d-9
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
      m = 30
c...     for cgstab(l):
      l= 2
c...     maximum number of matvecs in cgstab or gmres
      maxnmv = 100
      testspace = 3
c...     Testspace 1: w = "Standard Petrov" * v (Section 3.1.1)
c...     Testspace 2: w = "Standard 'variable' Petrov" * v (Section 3.1.2)
c...     Testspace 3: w = "Harmonic Petrov" * v (Section 3.5.1)
c
      if ( method .eq. 1 ) then
         lwork =  4 +  m  + 5*jmax + 3*kmax
      else
         lwork = 10 + 6*l + 5*jmax + 3*kmax
      end if
      wanted = .true.
c
      call jdqz(alpha, beta, eivec, wanted, n, target, tol, 
     $     kmax, jmax, jmin,
     $     method, m, l, maxnmv, maxstep,
     $     lock, order, testspace, zwork, lwork )
c
      elapse = etime( tarray )
c
c...     Compute the norms of the residuals:
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
      write(*,10) tarray(1), elapse
c
   10 format(1x,'END JDQZ AFTER ',f6.2,' SEC. CPU-TIME AND ', f6.2,
     $       ' SEC. ELAPSED TIME' )
      end

      subroutine PRECON( neq, q )
c...............................................
c...     Subroutine to compute q = K^-1 q
c...............................................
      integer        neq, i
      double complex q(neq)
c
      do i = 1, neq
         q(i) = i*q(i)/(i*i-31)
      end do
c
      end
c
      subroutine AMUL( neq, q, r )
c...............................................
c...     Subroutine to compute r = Aq
c...............................................
      integer        neq, i
      double complex q(neq), r(neq)
c
      do i = 1, neq
	 r(i) = i*q(i)
      end do
c
      end
c
      subroutine BMUL( neq, q, r )
c...............................................
c...     Subroutine to compute r = Bq
c...............................................
      integer        neq, i
      double complex q(neq), r(neq)
c
      do i = 1, neq
	 r(i) = q(i)/i
      end do
c
      end

