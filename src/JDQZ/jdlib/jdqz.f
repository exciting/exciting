      subroutine jdqz (alpha, beta, eivec, wanted,
     $     n, shift, eps, kmax, jmax, jmin,
     $     method, m, l,
     $     maxnmv, maxstep, lock, order, testspace, work, lwork)
c
c     Programmer: Diederik R. Fokkema
c     Modified         : M. van Gijzen
c     Modified 05-24-96: M. Kooper: ra and rb, the Schur matrices of A and B, 
c              added, as well as the vectors sa and sb, which contain the
c              innerproducts of ra with Z and rb with Z. This is added to be
c              enable to compute the eigenvectors in EIVEC
c     Modification 08-27-96: Different computation of eigenvectors, MvG
c
c     .. Parameters ..
c
      implicit none
      integer gmres, cgstab
      parameter ( gmres = 1, cgstab = 2 )
      integer kmax, jmax, jmin, method, m, l, maxnmv, maxstep, order
      integer testspace
      double precision eps, lock
      double complex shift
      integer n, lwork
      double complex work(n,*), eivec(n,*), alpha(*), beta(*)
      logical wanted
c
c     .. Local Parameters ..
c
      logical loop, try, found, ok
      integer ldvs, ldzwork, iseed(4)
      parameter (ldvs=50, ldzwork=4*ldvs)
      double complex ma(ldvs,ldvs), mb(ldvs,ldvs),
     $     zma(ldvs,ldvs), zmb(ldvs,ldvs),
     $     vsl(ldvs,ldvs), vsr(ldvs,ldvs),
     $     ra(ldvs,ldvs), rb(ldvs, ldvs),
     $     zwork(4*ldvs), aconv(ldvs), bconv(ldvs)

      integer ldqkz
      parameter (ldqkz=ldvs)
      integer ipivqkz(ldqkz)
      double complex mqkz(ldqkz,ldqkz), invqkz(ldqkz,ldqkz), f(ldqkz)

      integer i, j, k, info, mxmv, step
      integer d, u, v, w, aux, av, bv, q, z, kz, itmp, tp
      integer solvestep
      double precision rnrm, rwork(3*ldvs), deps
      double precision dtmp
      double complex zalpha, zbeta, targeta, targetb, evcond
      double complex shifta, shiftb
      double complex zero, one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
c
c     .. Functions ..
c
      double precision dznrm2
c
c     .. Data ..
c
      data iseed /3,3,1966,29/

c     ... Executable Statements
c
c...     Are there errors in the parameters?
      if ( kmax .gt. jmax )
     $   call error('jdqz: kmax greater than jmax')
      if ( jmax .gt. 50 )
     $   call error('jdqz: jmax greater than 50')
      if ( method .ne. 1 .and. method .ne. 2 )
     $   call error('jdqz: illegal choice for solution method')
      if ( order .lt. -2 .or. order .gt. 2 )
     $   call error('illegal value for order, must be between -2 and 2')

c...     d   = rhs, these pointers refer to the columns of the workspace
      d   = 1
c...     Workspace for jdqzmv
      tp  = d+1
c...     u   = pointer to Krylov space GMRES(m) or Bi-CSTAB(l)
      u   = tp+1
c...     v   = pointer to search space JDQZ with max dimension jmax
      if ( method .eq. gmres ) then
	 v   = u+m+1
      else if ( method .eq. cgstab ) then
	 v   = u+2*l+6
      end if
c...     w   = pointer to test subspace JDQZ with max dimension jmax
      w   = v+jmax
c...     av  = pointer to subspace AV with max dimension jmax
      av  = w+jmax
c...     bv  = pointer to subspace BV with max dimension jmax
      bv  = av+jmax
c...     aux =
      aux = bv+jmax
c...     q   = pointer to search Schur basis in JDQZ with max dimension kmax
      q   = aux+jmax
c...     z   = pointer to test Schur basis in JDQZ with max dimension kmax
      z   = q+kmax
c...     kz  = pointer to matrix K^{-1}Z_k
      kz  = z+kmax
      if (kz+kmax-1.gt.lwork) call error ('qz: memory fault')
c
c     --- initialize loop
c

      ok = .true.

      evcond = dcmplx(sqrt(abs(shift)**2+abs(one)**2))
      shifta = shift/evcond
      shiftb = one/evcond

      targeta = shifta
      targetb = shiftb

      zalpha = shifta
      zbeta = shiftb

      step = 0
      deps = dble(one)
      mxmv = 0

      solvestep = 0

      j = 0
      k = 0

c
c     --- loop
c
 100  continue
      loop = (k.lt.kmax.and.step.lt.maxstep)
      if (loop) then
	 step = step+1
	 solvestep = solvestep+1
	 if (j.eq.0) then
	    call zlarnv(2, iseed, n, work(1,v+j))
	    call zlarnv(2, iseed, n, work(1,w+j))
	    do i=1,n
	       dtmp = dble(work(i,v+j))
	       work(i,v+j) = dcmplx(dtmp,0d0)
	       dtmp = dble(work(i,w+j))
	       work(i,w+j) = dcmplx(dtmp,0d0)
	    enddo
	 else
	    mxmv = maxnmv
	    deps = 2d0**(-solvestep)
	    if (j.lt.jmin) then
	       mxmv = 1
	       call zgmres (n, work(1,v+j), work(1,d), m, deps,
     $              mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
	    elseif (method.eq.gmres) then
	       mxmv = m
	       call zgmres (n, work(1,v+j), work(1,d), m, deps,
     $              mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
	    elseif (method.eq.cgstab) then
	       call zcgstabl (n, work(1,v+j), work(1,d), l,
     $              deps, mxmv, zalpha,
     $              zbeta, k+1, work(1,kz), work(1,q), invqkz,
     $              ldqkz, ipivqkz, f, work(1,u), 2*l+6)
	    endif
	 endif
	 j = j+1

	 call zmgs (n, j-1, work(1,v), work(1,v+j-1), 1 )
	 call zmgs (n, k, work(1,q), work(1,v+j-1), 1 )

	 if (testspace.eq.1) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp), 
     $           -dconjg(shiftb), dconjg(shifta))
	 elseif (testspace.eq.2) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp),
     $           -dconjg(zbeta), dconjg(zalpha))
	 elseif (testspace.eq.3) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp),
     $           shifta, shiftb)
	 endif

         call zmgs (n, j-1, work(1,w), work(1,w+j-1), 1 )
 	 call zmgs (n, k, work(1,z), work(1,w+j-1), 1 )

         call amul(n, work(1,v+j-1), work(1,av+j-1))
         call bmul(n, work(1,v+j-1), work(1,bv+j-1))

         call makemm (n, j, work(1,w), work(1,av), ma, zma, ldvs)
         call makemm (n, j, work(1,w), work(1,bv), mb, zmb, ldvs)

         call zgegs ('v', 'v', j, zma, ldvs, zmb, ldvs,
     $        alpha, beta, vsl, ldvs, vsr, ldvs, zwork,
     $        ldzwork, rwork, info)

         try = .true.
 200     continue
         if (try) then
c
c           --- Sort the Petrov pairs ---
c
            call qzsort (targeta, targetb, j, zma, zmb, vsl, vsr,
     $           ldvs, alpha, beta, order)

            zalpha = zma(1,1)
            zbeta = zmb(1,1)

            evcond = dcmplx(sqrt(abs(zalpha)**2+abs(zbeta)**2))
c
c            --- compute new q ---
c
            call zgemv ('n', n, j, one, work(1,v), n, vsr(1,1),
     $           1, zero, work(1,q+k), 1)
            call zmgs (n, k, work(1,q), work(1,q+k), 1 )
c
c           --- compute new z ---
c
            call zgemv ('n', n, j, one, work(1,w), n, vsl(1,1),
     $           1, zero, work(1,z+k), 1)
            call zmgs (n, k, work(1,z), work(1,z+k), 1 )
c
c           --- Make new qkz ---
c
            call zcopy (n, work(1,z+k), 1, work(1,kz+k), 1)
            call precon (n, work(1,kz+k))
            call mkqkz (n, k+1, work(1,q), work(1,kz), mqkz, invqkz,
     $           ldqkz, ipivqkz)
c
c           --- compute new (right) residual= beta Aq - alpha Bq and
c               orthogonalize this vector on Z.
c
            call jdqzmv (n, work(1,q+k), work(1,d), work(1,tp),
     $           zalpha, zbeta)
            call zmgs (n, k, work(1,z), work(1,d), 0 )

            rnrm = dznrm2 (n, work(1,d), 1)/dble(evcond)

            if (rnrm.lt.lock.and.ok) then
               targeta = zalpha
               targetb = zbeta
               ok = .false.
            endif

            found = (rnrm.lt.eps.and.
     $           (j.gt.1.or.k.eq.kmax-1))
            try =  found

            if (found) then

c              --- increase the number of found evs by 1 ---
               k = k+1

c              --- store the eigenvalue
               aconv(k) = zalpha
               bconv(k) = zbeta

               solvestep = 0
               if (k.eq.kmax) goto 100
               call zgemm ('n', 'n', n, j-1, j, one, work(1,v), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = v
               v = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,av), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = av
               av = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,bv), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = bv
               bv = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,w), n,
     $              vsl(1,2), ldvs, zero, work(1,aux), n)
               itmp = w
               w = aux
               aux = itmp
               j = j-1
               call zlacpy ('a', j, j, zma(2,2), ldvs, ma, ldvs)
               call zlacpy ('a', j, j, ma, ldvs, zma, ldvs)
               call zlacpy ('a', j, j, zmb(2,2), ldvs, mb, ldvs)
               call zlacpy ('a', j, j, mb, ldvs, zmb, ldvs)
               call zlaset ('a', j, j, zero, one, vsr, ldvs)
               call zlaset ('a', j, j, zero, one, vsl, ldvs)
               targeta = shifta
               targetb = shiftb
               ok = .true.
               mxmv = 0
               deps = dble(one)
            else if (j.eq.jmax) then
               call zgemm ('n', 'n', n, jmin, j, one, work(1,v), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = v
               v = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,av), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = av
               av = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,bv), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = bv
               bv = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,w), n,
     $              vsl, ldvs, zero, work(1,aux), n)
               itmp = w
               w = aux
               aux = itmp
               j = jmin
               call zlacpy ('a', j, j, zma, ldvs, ma, ldvs)
               call zlacpy ('a', j, j, zmb, ldvs, mb, ldvs)
               call zlaset ('a', j, j, zero, one, vsr, ldvs)
               call zlaset ('a', j, j, zero, one, vsl, ldvs)
            endif
            goto 200
         endif
         goto 100
      endif

c
c...     Did enough eigenpairs converge?
      kmax = k

      if ( wanted ) then
c...        Compute the Schur matrices if the eigenvectors are
c...        wanted, work(1,tp) is used for temporary storage
c...        Compute RA:
         call zlaset ('l', k, k, zero, zero, ra, ldvs)    
         do i = 1, k
            call amul( n, work(1,q+i-1), work(1,tp) )
            call zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, ra(1,i), 1 )
         end do
c...        Compute RB:
         call zlaset ('l', k, k, zero, zero, rb, ldvs)    
         do i = 1, k
            call bmul( n, work(1,q+i-1), work(1,tp) )
            call zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, rb(1,i), 1 )
         end do

c        --- The eigenvectors RA and RB  belonging to the found eigenvalues
c            are computed. The Schur vectors in VR and VS are replaced by the
c            eigenvectors of RA and RB
         call zgegv('N','V',k,ra, ldvs, rb, ldvs, alpha, beta,vsl, ldvs,
     $              vsr,ldvs, zwork,ldzwork,rwork,info)
c        --- Compute the eigenvectors belonging to the found eigenvalues
c            of A and put them in EIVEC
         call zgemm('n', 'n', n, k, k, one, work(1,q), n,
     $              vsr, ldvs, zero, eivec, n)
      else
c
c...        Store the Schurvectors in eivec:
         call zcopy( k*n, work(1,q), 1, eivec, 1 )
         call zcopy( k, aconv, 1, alpha, 1 )
         call zcopy( k, bconv, 1, beta, 1 )
      end if

      end
