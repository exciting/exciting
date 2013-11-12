
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine iterativearpacksecequn (ik, ispn, apwalm, vgpc, evalfv, &
& evecfv)
!
      Use modinput
  !USES:
      Use modfvsystem
      Use modmpi
      Use mod_eigensystem
      Use mod_timing
      Use mod_Gkvector
      Use mod_potential_and_density
      Use mod_muffin_tin
      Use mod_atoms, Only: natmtot
!
      Use mod_spin, Only: nspnfv
      Use mod_APW_LO, Only: apwordmax
      Use mod_eigenvalue_occupancy, Only: nstfv
      Use sclcontroll
!
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   ispn   : first-variational spin index (in,integer)
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   vgpc   : G+k-vectors in Cartesian coordinates
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  ! This routine will perform several ARPACK iterations
!
  !BOC
      Implicit None
#ifdef DEBUG
  !include declarations for timing output of ARPACK
#include "./debugf90.h"
#endif
  ! arguments
      Integer, Intent (In) :: ik
      Integer, Intent (In) :: ispn
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Inout) :: evalfv (nstfv, nspnfv)
      Complex (8), Intent (Inout) :: evecfv (nmatmax, nstfv, nspnfv)
!
  ! local variables
      Logical :: packed
      Type (evsystem) :: system
!
      Integer :: n
      Real(8) :: cpu0, cpu1, cpu2
      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
  !IO vars
      Integer :: koffset, recl
      Character (256) :: outfilenamestring, filetag
      External outfilenamestring
!
  !ARPACK Interface vars
      Integer :: ido, nev, ncv, lworkl, info, infoznaupd, info2, j, i
      Integer :: nevmax, ncvmax, nmax
      Integer :: nconv, maxitr, ishfts, mode, ldv
      Integer :: iparam (11), ipntr (14)
      Complex (8), Allocatable :: resid (:), v (:, :), workev (:), &
     & workl (:), d (:)
      Complex (8), Pointer :: workd (:)
      Real (8), Allocatable :: rwork (:), rd (:)
      Integer, Allocatable :: idx (:)
      Complex (8) :: sigma
      Character :: bmat * 1, which * 2
      Real (8) :: tol
      Logical :: rvec
      Logical :: select (nmat(ispn, ik))
      Complex (8), Pointer :: vin (:), vout (:)

      logical :: ImproveInverse
      type (HermitianMatrix) :: zm
      Complex (8), Allocatable :: bres(:),vecupd(:),bvec(:)
      integer :: ii,ib
      real(8) :: berr
      
      logical :: spLU
      Complex (4),allocatable :: cm(:,:)
      integer, allocatable :: luipiv(:)
      integer :: luinfo
      real(8) ta,tb

      real(8) :: DegeneracyThr
      integer :: blockstart,blocksize
      Complex (8), Allocatable :: h(:,:),blockH(:,:),blockS(:,:)
      Complex (8), Allocatable :: zm2(:,:)
      real(8), allocatable :: evals(:)
      logical :: SeparateDegenerates,sorted

ImproveInverse=.false.
spLU=.false.
SeparateDegenerates=.true.
DegeneracyThr=1d-7
!
#ifdef DEBUG
      ndigit = - 3
      logfil = 6
      mngets = 1
      mnaitr = 1
      mnapps = 1
      mnaupd = 1
      mnaup2 = 1
      mneigh = 1
      mneupd = 1
      Open (logfil, File="ARPACK.OUT", Action="WRITE")
#endif
!
!
  !##################
  !ARPACK parameters
  !##################
      nev = nstfv
      ncv = 2 * nev
      ncv = Min (2*nev, maxncv)
      ncv = Max (ncv, nev+2)
      nevmax = nev
      ncvmax = ncv
      nmax = nmatmax
      n = nmat (ispn, ik)
      ldv = n
      lworkl = 3 * ncvmax * ncvmax + 5 * ncvmax
      Allocate (workd(3*nmax))
      Allocate (resid(nmax))
      Allocate (v(ldv, ncvmax))
      Allocate (workev(2*ncvmax))
      Allocate (workl(lworkl))
      Allocate (d(ncvmax))
      Allocate (rwork(ncvmax))
      Allocate (rd(ncvmax), idx(ncvmax))
      bmat = 'G'
      which = 'LM'
      if (input%groundstate%solver%ArpackUserDefinedShift) then
        sigma=input%groundstate%solver%ArpackShift
      else
        If (lowesteval .Eq.-1.d0) Then
          Call minenergy (sigma)
        Else
          sigma = dcmplx (lowesteval, 0)
        End If
      endif
      resid (:) = 0.0
      tol = input%groundstate%solver%epsarpack
      ido = 0
      info = 0
      ishfts = 1
      maxitr = 400 * nstfv
      mode = 3
      iparam (1) = ishfts
      iparam (3) = maxitr
      iparam (7) = mode
  !################################
  !open file with previous residual
  !################################
      Inquire (IoLength=Recl) resid
      koffset = ik - firstk (procofk(ik)) + 1
      infoznaupd = 0
!
  !##################
  !setup hamiltonian#
  !##################
!
!
      If (associated(input%groundstate%solver)) Then
         packed = input%groundstate%solver%packedmatrixstorage
      Else
         packed = .True.
      End If
!
      Call newsystem (system, packed, n)
      h1on=(input%groundstate%ValenceRelativity.eq.'lkh')
      Call hamiltonandoverlapsetup (system, ngk(ispn, ik), apwalm, &
     & igkig(1, ispn, ik), vgpc)
!
!
      Call timesec(cpu0)
  !#######################################################################
  !calculate LU decomposition to be used in the reverse communication loop
  !#######################################################################
!
      if (SeparateDegenerates) then 
        allocate(h(n,n))
        h=system%hamilton%za
      endif
     

      Call HermitianMatrixAXPY (-sigma, system%overlap, system%hamilton)

if (improveinverse) then
      call newmatrix(zm,.false.,system%overlap%rank)
      zm%za=system%hamilton%za
      allocate(bres(n))
      allocate(vecupd(n))
      allocate(bvec(n))
endif
      if (spLU) then
        allocate(cm(n,n))
        allocate(luipiv(n)) 
        allocate(system%hamilton%ipiv(n))
        cm(:,:)=cmplx(system%hamilton%za(1:n,1:n))
        Call CGETRF (n, n, cm, n, luipiv, luinfo)
        system%hamilton%za(1:n,1:n)=dcmplx(cm(1:n,1:n))
        system%hamilton%ipiv=luipiv
        system%hamilton%ludecomposed = .True.
        deallocate(cm,luipiv)
      else
         Call HermitianMatrixLU (system%hamilton)
      endif
      Call timesec(cpu1)
  !################################################
  !# reverse comunication loop of arpack library: #
  !################################################
!
      Do i = 1, maxitr
       call timesec(ta)
         Call znaupd (ido, bmat, n, which, nev, tol, resid, ncv, v, &
        & ldv, iparam, ipntr, workd, workl, lworkl, rwork, infoznaupd)
       call timesec(tb)
!       write(*,*) 'znaupd',tb-ta
         vin => workd (ipntr(1) :ipntr(1)+n-1)
         vout => workd (ipntr(2) :ipntr(2)+n-1)
!
!         write(*,*) ido
         If (ido .Eq.-1 .Or. ido .Eq. 1) Then
!
            call timesec(ta) 
            Call Hermitianmatrixvector (system%overlap, one, vin, zero, vout)

if (improveinverse) then
            bvec=vout
endif

            Call Hermitianmatrixlinsolve (system%hamilton, vout)

if (improveinverse) then
          berr=1d0
          do while(berr.gt.1d-6)
            Call Hermitianmatrixvector (zm, one, vout, zero, vecupd)
            bres=bvec-vecupd

            berr=0d0
            do ib=1,n
              berr=berr+abs(bres(ib))**2
            enddo
            berr=sqrt(berr/dble(n))

            Call Hermitianmatrixlinsolve (system%hamilton, bres)
            vout=vout+bres
          enddo
!          write(*,*)
!          read(*,*)
endif
          call timesec(tb)
!          write(*,*) 'ido-1',tb-ta


         Else If (ido .Eq. 1) Then
            Call zcopy (n, workd(ipntr(3)), 1, vout, 1)
if (improveinverse) then
            bvec=vout
endif
            Call Hermitianmatrixlinsolve (system%hamilton, vout)
if (improveinverse) then
          berr=1d0
          do while(berr.gt.1d-6)
            Call Hermitianmatrixvector (zm, one, vout, zero, vecupd)
            bres=bvec-vecupd

            berr=0d0
            do ib=1,n
              berr=berr+abs(bres(ib))**2
            enddo
            berr=sqrt(berr/dble(n))

            Call Hermitianmatrixlinsolve (system%hamilton, bres)
            vout=vout+bres
          enddo
!          write(*,*)
!          read(*,*)
endif

         Else If (ido .Eq. 2) Then
            call timesec(ta)
            Call Hermitianmatrixvector (system%overlap, one, vin, zero, vout)
            call timesec(tb)
!            write(*,*) 'Hermitianmatrixvector',tb-ta
         Else
            Exit
         End If 
!         read(*,*)
         
      End Do
  !###############
  ! errorhandling
  !###############
      If (infoznaupd .Ne. 0) Then
         Print *, ' '
         Print *, ' Error with znaupd, info = ', infoznaupd
         Print *, ' Check the documentation of znaupd'
         Print *, ' '
         Stop
      Else
!
         If (i .Gt. maxitr) Then
            Write (*,*) "Error reached maximum iteration count in arpac&
           &k."
            Stop
         End If
     !########################
     !post processing of evec
     !########################
         rvec = .True.
         select = .True.
         Call zneupd (rvec, 'A', select, d, v, n, sigma, workev, bmat, &
        & n, which, nev, tol, resid, ncv, v, n, iparam, ipntr, workd, &
        & workl, lworkl, rwork, info2)
         If (info2 .Ne. 0) Then
            Print *, ' '
            Print *, ' Error with zneupd, info = ', info2
            Print *, ' Check the documentation of zneupd'
            Print *, ' '
            Write (*,*) "eval", d (1:nev)
            Write (*,*) "iter", i
            Stop
         End If
!
      End If
      Call timesec(cpu2)
      timefv = timefv + cpu2 - cpu0
#ifdef DEBUG
      Close (logfil)
#endif
      If (rank .Eq. 0) write (60,*) "k=", ik, "ARPACK iterations", i
      If (rank .Eq. 0) write (60,*) "matrixsize", n, "time LU", cpu1 - &
     & cpu0, "iterations", cpu2 - cpu1
      If (rank .Eq. 0) write (60,*) "minenergy (inversioncenter)", dble &
     & (sigma)
!
  !##########################
  !sort and copy eigenvectors
  !##########################
      rd = real (d)
      Call sortidx (nstfv, rd(:), idx(:))
      Do j = 1, nstfv
!         write(*,*) d(j)
         evecfv (1:n, j, ispn) = v (1:n, idx(j))
         evalfv (j, ispn) = rd (idx(j))
      End Do
!      write(*,*)
      Deallocate (workd, resid, v, workev, workl, d)
      Deallocate (rwork, rd, idx)
      if (ImproveInverse) then
        call deletematrix(zm)
        deallocate(bres,vecupd,bvec)
      endif

if (SeparateDegenerates) then
        blockstart=1
        do i=2,nstfv
          if (evalfv (i, ispn)-evalfv (i-1, ispn).gt.DegeneracyThr) then
            blocksize=i-blockstart 
            if (blocksize.ge.2) then
              allocate(zm2(n,blocksize))
              allocate(blockH(blocksize,blocksize))
              allocate(blockS(blocksize,blocksize))
              allocate(workd(4*blocksize))
              allocate(rwork(3*blocksize-2))

              call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         n, &          ! M ... rows of op( A ) = rows of C
                         blocksize, &           ! N ... cols of op( B ) = cols of C
                         n, &          ! K ... cols of op( A ) = rows of op( B )
                         one, &          ! alpha
                         h, &           ! A
                         n,&           ! LDA ... leading dimension of A
                         evecfv(:,blockstart:blockstart+blocksize-1,ispn), &           ! B
                         nmatmax, &          ! LDB ... leading dimension of B
                         zero, &          ! beta
                         zm2, &  ! C
                         n &      ! LDC ... leading dimension of C
                         )
              call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         blocksize, &          ! M ... rows of op( A ) = rows of C
                         blocksize, &           ! N ... cols of op( B ) = cols of C
                         n, &          ! K ... cols of op( A ) = rows of op( B )
                         one, &          ! alpha
                         evecfv(:,blockstart:blockstart+blocksize-1,ispn), &           ! A
                         nmatmax,&           ! LDA ... leading dimension of A
                         zm2, &           ! B
                         n, &          ! LDB ... leading dimension of B
                         zero, &          ! beta
                         blockH, &  ! C
                         blocksize &      ! LDC ... leading dimension of C
                         )

              call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         n, &          ! M ... rows of op( A ) = rows of C
                         blocksize, &           ! N ... cols of op( B ) = cols of C
                         n, &          ! K ... cols of op( A ) = rows of op( B )
                         one, &          ! alpha
                         system%overlap%za, &           ! A
                         n,&           ! LDA ... leading dimension of A
                         evecfv(:,blockstart:blockstart+blocksize-1,ispn), &           ! B
                         nmatmax, &          ! LDB ... leading dimension of B
                         zero, &          ! beta
                         zm2, &  ! C       
                         n &      ! LDC ... leading dimension of C
                         )
              call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         blocksize, &          ! M ... rows of op( A ) = rows of C
                         blocksize, &           ! N ... cols of op( B ) = cols of C
                         n, &          ! K ... cols of op( A ) = rows of op( B )
                         one, &          ! alpha
                         evecfv(:,blockstart:blockstart+blocksize-1,ispn), &           ! A
                         nmatmax,&           ! LDA ... leading dimension of A
                         zm2, &           ! B
                         n, &          ! LDB ... leading dimension of B
                         zero, &          ! beta
                         blockS, &  ! C
                         blocksize &      ! LDC ... leading dimension of C
                         )

             call zhegv	(1, &                ! A*x = (lambda)*B*x
                         'V', &              ! Compute eigenvalues and eigenvectors
                         'U', &              ! Upper triangle
                         blocksize, &        ! Size of the problem
                         blockH, &           ! A
                         blocksize, &        ! leading dimension of A
                         blockS, &           ! B
                         blocksize, &        ! leading dimension of B
                         evalfv (blockstart:i-1, ispn), &                ! (lambda)
                         workd, &             ! work array
                         4*blocksize, &            ! size of lwork
                         rwork, &            ! another work array
                         info &              ! info
                        ) 	

              call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         n, &          ! M ... rows of op( A ) = rows of C
                         blocksize, &           ! N ... cols of op( B ) = cols of C
                         blocksize, &          ! K ... cols of op( A ) = rows of op( B )
                         one, &          ! alpha
                         evecfv(:,blockstart:blockstart+blocksize-1,ispn), &           ! B
                         nmatmax, &          ! LDB ... leading dimension of B
                         blockH, &           ! A
                         blocksize,&           ! LDA ... leading dimension of A
                         zero, &          ! beta
                         zm2, &  ! C
                         n &      ! LDC ... leading dimension of C
                         )
              evecfv(1:n,blockstart:blockstart+blocksize-1,ispn)=zm2(:,:)

              deallocate(zm2,blockS,blockH,workd,rwork)
            endif
            blockstart=i
          endif
        enddo
      deallocate(h)

endif
      
      Call deleteystem (system)


! We may have ruined the ascending order in the list of eigenvals and eigenvecs
! Restoring it if necessary
if (SeparateDegenerates) then
      sorted=.true.
      do i=2,nstfv
        sorted=(sorted.and.(evalfv(i,ispn).ge.evalfv(i-1,ispn)))
      enddo 
      if (.not.sorted) then  ! we sort the eigenvals and eigenvecs
        allocate(idx(nstfv))
        Call sortidx (nstfv, evalfv(:,ispn), idx(:))
        allocate(v(n,nstfv))        
        allocate(evals(nstfv))
        v(1:n,:)=evecfv(1:n,:,ispn)
        evals(:)=evalfv(:,ispn)
        Do j = 1, nstfv
          evecfv (1:n, j, ispn) = v (1:n, idx(j))
          evalfv (j, ispn) = evals (idx(j))
        End Do
        deallocate(v,evals,idx)
      endif
endif


! Nearly a replica of a part from seceqnfv
if (input%groundstate%ValenceRelativity.eq.'lkh') then
      Call newsystem (system, packed, n)
      h1aa=0d0
      h1loa=0d0
      h1lolo=0d0
      h1on=.false.
      Call hamiltonandoverlapsetup (system, ngk(ispn, ik), apwalm, &
     & igkig(1, ispn, ik), vgpc)
      call olprad
      allocate(h(n,nstfv))
      allocate(zm2(nstfv,nstfv))


      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n, &          ! M ... rows of op( A ) = rows of C
                  nstfv, &           ! N ... cols of op( B ) = cols of C
                  n, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  system%overlap%za, &           ! A
                  n,&           ! LDA ... leading dimension of A
                  evecfv(:,:,ispn), &           ! B
                  nmatmax, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  h, &  ! C
                  n &      ! LDC ... leading dimension of C
                  )
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  nstfv, &          ! M ... rows of op( A ) = rows of C
                  nstfv, &           ! N ... cols of op( B ) = cols of C
                  n, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  evecfv(:,:,ispn), &           ! A
                  nmatmax,&           ! LDA ... leading dimension of A
                  h, &           ! B
                  n, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  zm2, &  ! C
                  nstfv &      ! LDC ... leading dimension of C
                  )

      do i=1,nstfv
        evecfv(:,i,ispn)=evecfv(:,i,ispn)/sqrt(abs(zm2(i,i)))
      enddo
      deallocate(h,zm2)
      Call deleteystem (system)
endif
      Return
End Subroutine iterativearpacksecequn
!EOC
