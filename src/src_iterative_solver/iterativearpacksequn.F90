
! Copyright (C) 2005-2014 exciting team
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
      Use mod_kpoint, Only: nkpt
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
      type (HermitianMatrix) :: zm,cc
      Complex (8), Allocatable :: bres(:),vecupd(:),bvec(:)
      integer :: ii,ib,nounter
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

      real(8), allocatable :: rdiag(:)
      complex(8), allocatable :: cdiag(:),cinvdiag(:),zvec(:),rrvec(:),pvec(:),qvec(:)
      Complex (8), Allocatable :: work (:),cwork(:)
      integer, allocatable :: iwork(:)

      logical :: newarpackseed

      integer :: InvertMethod,counter
      character(len=4) :: strInvertMethod

      logical :: gev,keepworking
      integer :: method

      Type (HermitianMatrix), pointer ::pinverse
! sandbox variables :)
      complex(8), allocatable :: h3(:,:),zfft(:),h2(:,:)
      integer :: ig
      real(8) :: summa
      complex(8) :: rho,beta,alpha,rhoold
      logical :: usecg
 
      integer :: lwork
      complex(8) :: worksize
       

!     usecg=.true.
      gev=.true.
      SeparateDegenerates=.true.
      DegeneracyThr=1d-7
      spLU=(input%groundstate%solver%DecompPrec.eq.'sp')
      ImproveInverse=input%groundstate%solver%ArpackImproveInverse
      select case (input%groundstate%solver%ArpackLinSolve) 
      case ('LU')
        InvertMethod=LUdecomp
      case ('LDL')
        InvertMethod=LDLdecomp
      case ('LL')
        InvertMethod=LLdecomp
      case ('Diag')
        InvertMethod=Diagdecomp
      case ('InvertOnce')
        InvertMethod=InvertOnce
      end select
!     InvertMethod=InverseDiag
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
      info=0
      if (allocated(arpackseed)) then
        resid(1:n)=arpackseed(1:n,ik)
        infoznaupd=1
        arpackseed(:,ik)=0d0
        newarpackseed=.false.
      else
        resid (:) = 0.0
        infoznaupd = 0
        allocate(arpackseed(nmatmax,nkpt))
        arpackseed(:,:)=0d0
        newarpackseed=.true.
      endif
      tol = input%groundstate%solver%epsarpack
      ido = 0
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
      packed=.false.
!    
      Call newsystem (system, packed, n)
      h1on=(input%groundstate%ValenceRelativity.eq.'iora*')
      Call hamiltonandoverlapsetup (system, ngk(ispn, ik), apwalm, igkig(1, ispn, ik), vgpc)
!      write(*,*) '************ARPACK**************'
!      read(*,*)

! sandbox starts here
if (.false.) then
      Call HermitianMatrixAXPY (-sigma, system%overlap, system%hamilton)
      do i=1,100 !ngk(ispn,ik) 
!       write(*,*) abs(system%hamilton%za(i,i))
       summa=0d0
       do j=1,ngk(ispn,ik)
         summa=summa+abs(system%hamilton%za(j,i))**2
!         write(*,*) abs(system%hamilton%za(j,i))
       enddo
       write(*,*) abs(system%hamilton%za(i,i))/sqrt(summa)
      enddo
 
      stop
      write(*,*) ngkfft
      allocate(h2(ngktotfft,ngktotfft))
      allocate(h3(ngktotfft,ngk(ispn,ik)))
      allocate(zfft(ngktotfft))
      h2=0d0
      h3=0d0
      do j=1,ngk(ispn,ik)
       zfft=0d0
       do i=1,ngk(ispn,ik)
        ig=igkfft(i,ik)
        zfft(ig)=system%hamilton%za(i,j)
       enddo
       Call zfftifc (3, ngkfft, 1, zfft)
       h3(:,j)=zfft(:)
      enddo
      do j=1,ngktotfft
       zfft=0d0
       do i=1,ngk(ispn,ik)
        ig=igkfft(i,ik)
        zfft(ig)=h3(j,i)
       enddo
       Call zfftifc (3, ngkfft, -1, zfft)
       h2(j,:)=zfft(:)
      enddo
      write(*,*)
      do i=1,10!ngkfft(1)/2,ngkfft(1)/2 !ngktotfft/2,ngktotfft/2 !ngktotfft
       summa=0d0
       do j=1,ngktotfft
         summa=summa+abs(h2(j,i))**2
!        write(*,*) abs(h2(j,i)) !sqrt(abs(h3(j,i)**2)/(abs(h3(i,i)*abs(h3(j,j)))))
       enddo
       write(*,*) abs(h2(i,i)),sqrt(summa)
      enddo
      stop
!      Call zfftifc (3, ngkfft, 1, zfft)
endif
! end of sandbox

if (gev) then
       
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
      if (ImproveInverse) then
        zm%sp=.false.
        call newmatrix(zm,.false.,n)
        zm%za=system%hamilton%za
        allocate(bres(n))
        allocate(vecupd(n))
        allocate(bvec(n))

        allocate(pvec(n))
        allocate(zvec(n))
        allocate(qvec(n))
        allocate(rrvec(n))
        rho=0d0
      endif

if (InvertMethod.eq.InverseDiag) then
      allocate(cdiag(n))
      allocate(cinvdiag(n))
      do i=1,n
       cdiag(i)=system%hamilton%za(i,i)
       cinvdiag(i)=1d0/system%hamilton%za(i,i)
      enddo
      allocate(pvec(n))
      allocate(zvec(n))
      allocate(qvec(n))
      allocate(rrvec(n))
elseif (InvertMethod.eq.InvertOnce) then
      if (.not.associated(arpackinverse)) then
        write(*,*) 'init arpackinverse'
        allocate(arpackinverse(nkpt)) 
        do i=1,nkpt
          nullify(arpackinverse(i)%za)
        enddo
      endif
      
      pinverse=>arpackinverse(ik)
      if (.not.associated(pinverse%za)) then
        pinverse%sp=.false.
!        if (.not.associated(pinverse%za)) then
          call newmatrix(pinverse,.false.,n)
!        else
!          write(*,*) 'sum',sum(pinverse%za)
!          call deletematrix(pinverse)
!          call newmatrix(pinverse,.false.,n)
!        endif
!        write(*,*) pinverse%rank,n
        pinverse%za(1:n,1:n)=system%hamilton%za(1:n,1:n)
        method=LLDecomp
!        pinverse%ludecomposed=.false.
        call HermitianMatrixFactorize(pinverse,method)
        call HermitianMatrixInvert(pinverse,method)
      endif
!         cc%sp=.false.
!         call newmatrix(cc,.false.,n)
!         cc%za=0d0
!         call HermitianMatrixMatrix(cc,pinverse%za,system%hamilton%za,n,n,n)
!         write(*,*) cc%za(1:3,1:3)
!         call deletematrix(cc)
!         write(*,*) 'sum',sum(pinverse%za)
!      read(*,*) 
      system%hamilton%za=pinverse%za
      system%hamilton%ludecomposed=.true.
      
else
      if (InvertMethod.ne.Diagdecomp) then
        if (spLU) then
          cc%sp=.true.
          call newmatrix(cc,.false.,n)
          cc%ca=cmplx(system%hamilton%za)
          call HermitianMatrixFactorize(cc,InvertMethod)
          system%hamilton%za=dcmplx(cc%ca)
          if (associated(cc%ipiv)) then
            allocate(system%hamilton%ipiv(n))
            system%hamilton%ipiv=cc%ipiv
          endif
          system%hamilton%ludecomposed=.true.
          call deletematrix(cc)
        else
          Call HermitianMatrixFactorize(system%hamilton,InvertMethod)
        endif
      else
!       system%hamilton%za=system%overlap%za
        allocate(rdiag(n))
        allocate(cdiag(n))
        allocate(work(2*n-1))
        allocate(rwork(3*n-2))
        call zheev ('V', &                     ! compute eigenvals and eigenvecs
                    'U', &                     ! upper triangle is stored
                     n, &                       ! size of the matrix A
!                     system%overlap%za, &      ! A
                     system%hamilton%za, &      ! A
                     n, &                       ! leading dimension
                     rdiag, &                   ! eigenvals
                     work, &                    ! work array
                     2*n-1, &                 ! size of work
                     RWORK, &                    ! another work array
                     info &                     ! error message
                     )
        write(*,*) rdiag
        stop
        deallocate(work,rwork)
        

        do i=1,n
          cdiag(i)=dcmplx(1d0/rdiag(i),0d0)
        enddo
        deallocate(rdiag)
        allocate (zvec(n))

      endif
endif
      
      Allocate (rwork(ncvmax))

      Call timesec(cpu1)
!      write(*,*) 'LU',cpu1-cpu0
  !################################################
  !# reverse comunication loop of arpack library: #
  !################################################
!
      Do i = 1, maxitr
       call timesec(ta)
         Call znaupd (ido, bmat, n, which, nev, tol, resid, ncv, v, &
        & ldv, iparam, ipntr, workd, workl, lworkl, rwork, infoznaupd)
       call timesec(tb)
      
         vin => workd (ipntr(1) :ipntr(1)+n-1)
         vout => workd (ipntr(2) :ipntr(2)+n-1)

         If (ido .Eq.-1 .Or. ido .Eq. 1) Then
!***********************
! Linear solve segment *
! (H-sigma*S)*vout=vin *
!***********************
            call timesec(ta) 

            Call Hermitianmatrixvector (system%overlap, one, vin, zero, vout)
            if (ImproveInverse) then
              bvec=vout
            endif

If (InvertMethod.ne.InverseDiag) then
            if (InvertMethod.ne.Diagdecomp) then
              
              Call Hermitianmatrixlinsolve (system%hamilton, vout,InvertMethod)

            else
              Call zgemv ('C', n, n, one, system%hamilton%za, n, vout, 1, zero, zvec, 1)
              zvec(:)=zvec(:)*cdiag(:)
              Call zgemv ('N', n, n, one, system%hamilton%za, n, zvec, 1, zero, vout, 1)
            end if
else 
              rrvec=vout 
              zvec(:)=cinvdiag(:)*rrvec(:)
              rho=0d0
              do ib=1,n
                rho=rho+conjg(rrvec(ib))*zvec(ib)
              enddo
              pvec=zvec
              Call Hermitianmatrixvector (system%hamilton, one, pvec, zero, qvec)
              
              alpha=0d0
              do ib=1,n
                alpha=alpha+conjg(pvec(ib))*qvec(ib)
              enddo   
              alpha=rho/alpha
              vout=vout+alpha*pvec
              rrvec=rrvec-alpha*qvec
              bres=rrvec
            
endif

! Iterative improvement of linear solve if necessary
            if (ImproveInverse) then
If (invertMethod.ne.InverseDiag) then
!              bres=bvec
!              Call Hermitianmatrixvector (zm, -one,vout, one, bres)
                Call Hermitianmatrixvector (zm, one, vout, zero, vecupd)
                bres=bvec-vecupd
                rrvec=bres

!                rrho=0d0
!                do ib=1,n
!                  rho=rho+conjg(rrvec(ib))*zvec(ib)
!                enddo

endif
              berr=1d0
              counter=0
              do while ((berr.gt.1d-10).and.(counter.lt.10))
                berr=0d0
                do ib=1,n
                  berr=berr+abs(bres(ib))**2
                enddo
             
                berr=sqrt(berr/dble(n))
!                write(*,*) counter,berr,InvertMethod
!                read(*,*)
                if (berr.gt.1d-10) then                
!                  read(*,*)
                  counter=counter+1
If (invertMethod.ne.InverseDiag) then
                  zvec=rrvec
                  if (InvertMethod.ne.Diagdecomp) then
                    Call Hermitianmatrixlinsolve (system%hamilton, zvec,InvertMethod)
                  else
!                    Call zgemv ('C', n, n, one, system%hamilton%za, n, rrvec, 1, zero, zvec, 1)
!                    zvec(:)=zvec(:)*cdiag(:)
!                    Call zgemv ('N', n, n, one, system%hamilton%za, n, zvec, 1, one, vout, 1)
                  end if
                  rhoold=rho
                  rho=0d0
                  do ib=1,n
                    rho=rho+conjg(rrvec(ib))*zvec(ib)
                  enddo
                  if (counter.eq.1) then
                    pvec=zvec
                  else
                    beta=rho/rhoold
                    pvec=zvec+beta*pvec
                  endif
                  Call Hermitianmatrixvector (zm, one, pvec, zero, qvec)

                  alpha=0d0
                  do ib=1,n
                    alpha=alpha+conjg(pvec(ib))*qvec(ib)
                  enddo
                  alpha=rho/alpha
                  vout=vout+alpha*pvec
                  rrvec=rrvec-alpha*qvec
                  bres=rrvec

!                  if (InvertMethod.ne.Diagdecomp) then
!                    Call Hermitianmatrixlinsolve (system%hamilton, bres,InvertMethod)
!                    vout=vout+bres
!                  else
!                    Call zgemv ('C', n, n, one, system%hamilton%za, n, bres, 1, zero, zvec, 1)
!                    zvec(:)=zvec(:)*cdiag(:)
!                    Call zgemv ('N', n, n, one, system%hamilton%za, n, zvec, 1, one, vout, 1)
!                  end if
!                  Call Hermitianmatrixvector (zm, one, vout, zero, vecupd)
!                  bres=bvec-vecupd
else
                   zvec(:)=rrvec(:)*cinvdiag(:)
                   rhoold=rho
                   rho=0d0
                   do ib=1,n
                    rho=rho+conjg(rrvec(ib))*zvec(ib)
                   enddo
                   beta=rho/rhoold
                   pvec=zvec+beta*pvec
                   Call Hermitianmatrixvector (system%hamilton, one, pvec, zero, qvec)

                   alpha=0d0
                   do ib=1,n
                     alpha=alpha+conjg(pvec(ib))*qvec(ib)
                   enddo
                   alpha=rho/alpha
                   vout=vout+alpha*pvec
                   rrvec=rrvec-alpha*qvec 
                   bres=rrvec                  
endif
                endif
              enddo
              if ((counter.ge.10).and.(InvertMethod.eq.InvertOnce)) then
                pinverse=>arpackinverse(ik)
                pinverse%za=zm%za
                method=LLDecomp
                pinverse%ludecomposed=.false.
                call HermitianMatrixFactorize(pinverse,method)
                call HermitianMatrixInvert(pinverse,method)
                system%hamilton%za=pinverse%za
!                system%hamilton%ludecomposed=.true.
                vout=bvec
                Call Hermitianmatrixlinsolve (system%hamilton, vout,InvertMethod)
                Call Hermitianmatrixvector (zm, one, vout, zero, vecupd)
                bres=bvec-vecupd
                berr=0d0
                do ib=1,n
                  berr=berr+abs(bres(ib))**2
                enddo

                berr=sqrt(berr/dble(n))

                write(*,*) 'inverse updated' , berr
              endif

!              endif
!              write(*,*) 
            endif

            call timesec(tb)
!            write(*,*) tb-ta
!           read(*,*)

         Else If (ido .Eq. 2) Then
!******************************
! Matrix times vector segment *
! vout=S*vin                  *
!******************************
            call timesec(ta)
            Call Hermitianmatrixvector (system%overlap, one, vin, zero, vout)
            call timesec(tb)
         Else
            Exit
         End If 
         
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
!      write(*,*) 'iterations',cpu2-cpu1
      timefv = timefv + cpu2 - cpu0
#ifdef DEBUG
      Close (logfil)
#endif
      If (rank .Eq. 0) write (60,*) "k=", ik, "ARPACK iterations", i
      select case (InvertMethod) 
      case (LUdecomp)
        strInvertMethod="LU"
      case (LDLdecomp)
        strInvertMethod="LDL"
      case (LLdecomp)
        strInvertMethod="LL"
      case (Diagdecomp)
        strInvertMethod="Diag"
      case (InvertOnce)
        strInvertMethod="Inv"
      end select
    
     If (rank .Eq. 0) write (60,*) "matrixsize", n, "time ", strInvertMethod, cpu1 - &
     & cpu0, "iterations", cpu2 - cpu1
      If (rank .Eq. 0) write (60,*) "minenergy (inversioncenter)", dble &
     & (sigma)
!
  !##########################
  !sort and copy eigenvectors
  !##########################
 !     write(*,*) i
      rd = real (d)
      Call sortidx (nstfv, rd(:), idx(:))
      Do j = 1, nstfv
!         write(*,*) rd(idx(j))
         evecfv (1:n, j, ispn) = v (1:n, idx(j))
         evalfv (j, ispn) = rd (idx(j))
      End Do
      Deallocate (workd, resid, v, workev, workl, d)
      Deallocate (rwork, rd, idx)
     
      if (InvertMethod.eq.Diagdecomp) deallocate(zvec)
      if (ImproveInverse) then
        call deletematrix(zm)
        deallocate(bres,vecupd,bvec)
      endif
if (SeparateDegenerates) then
        blockstart=1
        do i=2,nstfv
          if (evalfv (i, ispn)-evalfv (i-1, ispn).gt.DegeneracyThr) then
            blocksize=i-blockstart 
            if (blocksize.ge.1) then
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
      Call deletesystem (system)
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
      Call deletesystem (system)
endif
!      do i=1,nstfv
!        arpackseed(1:n,ik)=arpackseed(1:n,ik)+evecfv(1:n,i,ispn)
!      enddo
      if (newarpackseed) then
        do i=1,nstfv
          arpackseed(1:n,1)=arpackseed(1:n,1)+evecfv(1:n,i,ispn)
        enddo
        do i=2,nkpt
          arpackseed(1:n,i)=arpackseed(1:n,1)
        enddo
      else
        do i=1,nstfv
          arpackseed(1:n,ik)=arpackseed(1:n,ik)+evecfv(1:n,i,ispn)
        enddo
        arpackseed(1:n,ik)=arpackseed(1:n,ik)/dble(nstfv)
      endif
! write(*,*) '**********ARPACK DONE*************'
else
!**************************************************************************
! This segment starts working if gev=.false.                              *
! Its purpose is to test what if we solve H*x=E*S*x by reducing it        *
! to the regular eigenvalue problem.                                      *
! Step 1. Decompose S=U**t *U , where t stands for dagger.                *
! Step 2. Calculate inv(U).                                               *
! Step 3. Find eigenvalues and eigenvectors for inv(U)**t *H*inv(U)*y=Ey  *
! Step 4 (not implemented). Calculate x=inv(U)y.                          *
! Step 5 (not implemented). Postprocess x - normalization, sorting etc.   *
!-------------------------------------------------------------------------*
! Observations.                                                           *
! a) Inverse of triangular matrix U costs about the same as the Cholesky  *
! decomposition. Very good!                                               *
! b) Iterations inv(U)**t *H*inv(U)*y cost 3 matrix-vector multiplies,    *
! which makes this approach comparable and typically more expensive       *
! time-wise to the basic shift-and-invert. Disappointing.                 *
! c) The number of iterations required is typically higher than in shift- *
! and-invert. Disappointing again.                                        *
! Conclusion. I keep the thing in the code just because it is curious,    *
! and I do not see a point in finishing it.                               *
! Cheers,                                                                 *
! Andris                                                                  *
!**************************************************************************

      Call HermitianMatrixAXPY (-sigma, system%overlap, system%hamilton)
      sigma=0d0
      allocate(zvec(n))
      allocate(bres(n))
      Allocate (rwork(ncvmax))
      
      ishfts = 1
      mode   = 1
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode 
!      system%hamilton%za=system%hamilton%za-sigma*system%overlap%za
!      Call HermitianMatrixAXPY (-sigma, system%overlap, system%hamilton)
if (.false.) then
      zm%sp=.false.
      call newmatrix(zm,.false.,n)
      zm%za=system%overlap%za
      call timesec(ta)
      call HermitianMatrixLL(zm)
      call timesec(tb)
      write(*,*) 'LL',tb-ta
      do i=1,n-1
         zm%za(i+1:n,i)=0d0
      enddo
      call timesec(ta)
      call ztrtri('U', &
                  'N', &
                   n,  &
                   zm%za, &
                   n, &
                   info & 
                  )        
      call timesec(tb)
endif
      write(*,*) 'inv',tb-ta       
      bmat='I'      
!     which='SR'
      which='LM'

      call timesec(ta)
      i=0
      keepworking=.true.
      do while ((i.lt.maxitr).and.keepworking)
        i=i+1
!        call timesec(ta)
        Call znaupd (ido, bmat, n, which, nev, tol, resid, ncv, v, &
          & ldv, iparam, ipntr, workd, workl, lworkl, rwork, infoznaupd)
!        call timesec(tb)
        vin => workd (ipntr(1) :ipntr(1)+n-1)
        vout => workd (ipntr(2) :ipntr(2)+n-1)

        If (ido .Eq.-1 .Or. ido .Eq. 1) Then
!            call timesec(ta)


!            Call zgemv ('N', n, n, one, zm%za, n, vin, 1, zero, zvec, 1)
!            Call zgemv ('N', n, n, one, system%hamilton%za, n, zvec, 1, zero, bres, 1)
!            Call zgemv ('C', n, n, one, zm%za, n, bres, 1, zero, vout, 1)
             Call zgemv ('N', n, n, one, system%hamilton%za, n, vin, 1, zero, vout, 1)

!            call timesec(tb)
        else 
          keepworking=.false.
        endif
      enddo
      call timesec(tb)
      write(*,*) 'iterations',tb-ta

      write(*,*) i,ido
      if (info.ne.0) then 
        write(*,*) 'ARPACK (znaupd) produced info=',info
        stop
      endif
      if ((ido.eq.-1).or.(ido.eq.1)) then
        write(*,*) 'Too many ARPACK iterations'
        stop
      endif

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
      write(*,*)
      do i=1,nstfv
        write(*,*) dble(d(i))
      enddo
      stop

! In case if someone decides to finish the solver,
! here is the spot for the postprocessing.
        
endif
      Return
End Subroutine iterativearpacksecequn
!EOC
