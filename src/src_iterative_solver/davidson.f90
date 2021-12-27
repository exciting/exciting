Subroutine davidson (system, nst, evecfv, evalfv,ik)
      Use modinput
  !USES:
      Use modfvsystem
      Use modmpi
      use modmain
      Use mod_eigensystem
      Use mod_timing
      Use mod_Gkvector
      Use mod_potential_and_density
      Use mod_muffin_tin
      Use mod_atoms, Only: natmtot
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
      Real (8), Intent (Inout) :: evalfv (nst)
      Complex (8), Intent (Inout) :: evecfv (nmatmax, nst)
      Type (evsystem) :: system
      integer, intent(in) :: nst
      integer, intent(inout) :: ik
!
  ! local variables
      Logical :: packed
      Type (HermitianMatrix) :: fm
!
      Integer :: n
      Real(8) :: cpu0, cpu1, cpu2, tsa,tsb
      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
! IO vars
      Integer :: koffset, recl
      Character (256) :: outfilenamestring, filetag
      External outfilenamestring
!
! ARPACK Interface vars
      Integer :: ido, nev, ncv, lworkl, info, infoznaupd, info2, j, i,it
      Integer :: nevmax, ncvmax, nmax, nstart
      Integer :: nconv, maxitr, ishfts, mode, ldv
      Integer :: iparam (11), ipntr (14)
      Complex (8), Allocatable :: resid (:), v (:, :), workev (:), &
     & workl (:), d (:)
      Complex (8), Pointer :: workd (:)
      Real (8), Allocatable :: rwork (:), rd (:), residlen(:)
      Integer, Allocatable :: idx (:)
      Complex (8) :: sigma
      Character :: bmat * 1, which * 2
      Real (8) :: tol,maxresid
      Logical :: rvec
      Logical :: select (nmat(1, ik))
      Complex (8), Pointer :: vin (:), vout (:)

      logical :: ImproveInverse
      type (HermitianMatrix) :: zm,cc
      Complex (8), Allocatable :: bres(:),vecupd(:)
      integer :: ii,ib,nounter
      real(8) :: berr
      
      real(8) :: DegeneracyThr
      integer :: blockstart,blocksize
      Complex (8), Allocatable :: h(:,:),blockH(:,:),blockS(:,:),Hx(:,:),Sx(:,:)
      Complex (8), Allocatable :: zm2(:,:),diag(:),zm3(:,:),zm4(:,:)
      Complex (8), allocatable :: trialvec(:,:),bufrecv(:)
      real(8), allocatable :: evals(:),sortedevals(:)
      logical :: SeparateDegenerates,sorted

      real(8), allocatable :: rdiag(:)
      complex(8), allocatable :: cdiag(:),cinvdiag(:),zvec(:),rrvec(:),pvec(:),qvec(:),extraS(:,:),extraH(:,:),extraSx(:,:),extraHx(:,:),ritzvec(:,:)
      Complex (8), Allocatable :: work (:),cwork(:), zax(:,:), sdiag(:),hdiag(:) 
      integer, allocatable :: iwork(:)

      logical :: newarpackseed

      integer :: counter,offset,gksize
      character(len=4) :: strInvertMethod

      logical :: gev,keepworking
      integer :: method

      Type (HermitianMatrix), pointer ::pinverse
! sandbox variables :)
      complex(8), allocatable :: h3(:,:),zfft(:),h2(:,:)
      integer :: ig
      real(8) :: summa,maxdiff,garums
      complex(8) :: rho,beta,alpha,rhoold, zsum,zsum2
      logical :: usecg, calcsingular
 
      integer :: lwork
      complex(8) :: worksize

      real(8) :: oldsum,newsum,norm,sft,newrd,maxkin,potl,time1,time2
      integer :: ndiv,nblocks,calls,npw,is,ia,blo,bhi,l,lm,ias,nsize,nadd,ilo,m,if1,nloall,n_local,npw_local,nusedsingular

      complex(8), external :: zdotc
      integer, allocatable :: gkoffset(:),ngklist(:),ibuf(:)

!call memory_report
write(*,*) 'ik=',ik
call timesec(tsa)

      tol = input%groundstate%solver%epsarpack
      packed=.false.

      npw=ngk(1,ik)
      n=nmat(1,ik)   !system%hamilton%rank
      nloall=n-npw

      npw_local=npw !gkhi_local-gklo_local+1
      n_local=nloall+npw_local

allocate(sdiag(n_local))
allocate(hdiag(n_local))
sdiag=0d0
hdiag=0d0

if (associated(system%hamilton%za)) then
do i=1,n
 hdiag(i)=system%hamilton%za(i,i)
 sdiag(i)=system%overlap%za(i,i)

enddo


else
sdiag(1:npw)= cfunig0

do i=1,npw
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      sdiag(i)=sdiag(i)+zdotc(mt_hscf%maxaa,system%apwi(1,i,ias),1,system%apwi(1,i,ias),1)
    enddo
  enddo
enddo

sdiag(npw_local+1:n_local)=1d0

hdiag(1:npw)= veffig0

if (input%groundstate%ValenceRelativity.ne."none") then
  do i=1,npw
    hdiag(i)=hdiag(i)+0.5d0*dot_product (current_vgkc(:, i), current_vgkc(:, i))*meffig(1)
  enddo
else
  do i=1,npw
    hdiag(i)=hdiag(i)+0.5d0*dot_product (current_vgkc(:, i), current_vgkc(:, i))*cfunig0
  enddo
endif
allocate(zvec(mt_hscf%maxaa))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(none) PRIVATE(i,is,ia,ias,zvec) SHARED(npw,nspecies,natoms,idxas,system,hdiag,mt_hscf)
!$OMP DO 
#endif
do i=1,npw
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      call zgemv('N',mt_hscf%maxaa,mt_hscf%maxaa,one,mt_hscf%main%aa(1,1,ias),mt_hscf%maxaa,system%apwi(1,i,ias),1,zero,zvec,1)
      hdiag(i)=hdiag(i)+zdotc(mt_hscf%maxaa,system%apwi(1,i,ias),1,zvec,1)
    enddo
  enddo
enddo

#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
deallocate(zvec)

 i=npw_local+1
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      if1=1
      do ilo=1,nlorb(is)
        l=lorbl(ilo, is)
        do m=-l,l
          hdiag(i+m+l)=mt_hscf%main%lolo(if1+m+l,if1+m+l,ias)
        enddo
        if1=if1+2*l+1
        i=i+2*l+1
      enddo
    enddo
  enddo

endif


      nusedsingular=nsingular 
      ndiv=nst
      nblocks=12
      calls=0
      allocate(rd(nblocks*ndiv+nusedsingular+nloall))
      allocate(trialvec(n_local,nblocks*ndiv+nusedsingular+nloall))
      trialvec=zzero

      allocate(ritzvec(n_local,ndiv))
      allocate(residlen(ndiv))


      allocate(Hx(n_local,nblocks*ndiv+nusedsingular+nloall))
      allocate(Sx(n_local,nblocks*ndiv+nusedsingular+nloall))
      allocate(extraS(nblocks*ndiv+nusedsingular+nloall,nblocks*ndiv+nusedsingular+nloall))
      allocate(extraH(nblocks*ndiv+nusedsingular+nloall,nblocks*ndiv+nusedsingular+nloall))
      allocate(extraSx(n_local,ndiv))
      allocate(extraHx(n_local,ndiv))


!--------------------------
! Subspace initialisation
!
       nsize=0
       nadd=0
Hx=zzero
Sx=zzero
! -> 1. all local orbitals
       do i=1,nloall
         nadd=nadd+1
         trialvec(i+npw_local,nsize+nadd)=1d0
       enddo
       nsize=nsize+nadd

       call HloSlo(n_local,npw_local,nsize,system,trialvec,Hx,Sx)
       if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
        call pacelo(n_local,npw_local,nsize,nstsv,0.25d0,pace(:,:,ik),trialvec ,Hx )
       endif

! -> 2. singular components
       if (nusedsingular.ne.0) then
         trialvec(1:npw_local,nsize+1:nsize+nusedsingular)=singular(1:npw_local,1:nusedsingular,ik)
         nsize=nsize+nusedsingular
         call GSortho(n_local,0,nsize-nloall,trialvec(:,nloall+1:))
       endif

! -> 3. initial guess for wavefunctions
       if (dble(sum(evecfv)).ne.0d0) then
!       wavefunctions from previous scf iterations
         trialvec(1:npw,nsize+1:nsize+nst)=evecfv(1:npw, 1:nst)
         nsize=nsize+nst
       else
!       or just a guess
         do i=1,ndiv
           trialvec(i,nsize+i)=one
           trialvec(i+1,nsize+i)=0.5d0
           trialvec(i+2,nsize+i)=0.25d0
         enddo
         nsize=nsize+ndiv
       endif

! GEV

call timesec(time1)

      allocate(BlockS(nsize,nsize))
      allocate(BlockH(nsize,nsize))
      BlockH=zero
      BlockS=zero
!
      call GSortho(n_local,nsize-ndiv-nloall,nsize-nloall,trialvec(:,nloall+1:))
      call HapwSapw(n_local,npw,nsize-nloall,system,trialvec(:,nloall+1:nsize),Hx(:,nloall+1:nsize),Sx(:,nloall+1:nsize),.true.)
      if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
write(*,*) 'pace activated'
        call paceapw(n_local,npw,nsize-nloall,nstsv,0.25d0,pace(:,:,ik),trialvec(1:n_local,nloall+1:nsize) ,Hx(1:n_local,nloall+1:nsize))
      endif

      call innerproduct(n_local,nsize,nsize,trialvec,Hx(:,1:nsize),blockH(:,1:nsize))
      call innerproduct(n_local,nsize,nsize,trialvec,Sx(:,1:nsize),blockS(:,1:nsize))
!      blockH(1:nsize-ndiv,1:nsize-ndiv)=extraH(1:nsize-ndiv,1:nsize-ndiv)
!      do i=1,nsize-ndiv
!        do j=nsize-ndiv+1,nsize
!          blockH(j,i)=conjg(blockH(i,j))
!        enddo
!      enddo
      extraH(1:nsize,1:nsize)=BlockH(1:nsize,1:nsize)

!      blockS(1:nsize-ndiv,1:nsize-ndiv)=extraS(1:nsize-ndiv,1:nsize-ndiv)
!      do i=1,nsize-ndiv
!        do j=nsize-ndiv+1,nsize
!          blockS(j,i)=conjg(blockS(i,j))
!        enddo
!      enddo
      extraS(1:nsize,1:nsize)=BlockS(1:nsize,1:nsize)

! The initial subspace is set up
!----------------------------------

!write(*,*)
!do i=1,nsize
!  write(*,*) blockH(i,i)
!enddo
!write(*,*) blockS(1,1)

      call diagH2(nsize,blockH,blockS,rd,info)
!if (info.eq.0) write(*,*) 'info=',info
!write(*,*)
!write(*,*) rd(1:ndiv)
!stop
      call pickritzvectors(nsize,ndiv,rd,nstart)
      call getritzvectors(n_local,nsize,ndiv,trialvec,blockH,ritzvec,nstart,npw_local)
!read(*,*)
#ifdef garums
      allocate(zm3(n,nsize))
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  nusedsingular, &          ! M ... rows of op( A ) = rows of C
                  nsize, &           ! N ... cols of op( B ) = cols of C
                  n, &          ! K ... cols of op( A ) = rows of op( B )
                  one, &          ! alpha
                  singular(1,1,ik), &           ! A
                  n,&           ! LDA ... leading dimension of A
                  zm3, &           ! B
                  n, &          ! LDB ... leading dimension of B
                  zero, &          ! beta
                  zm2, &  ! C
                  nusedsingular &      ! LDC ... leading dimension of C
                  )
      garums=sqrt(abs(zdotc(n,zm2(1,1),1,zm2(1,1),1)))
      do i=1,nusedsingular
        write(*,*) abs(zm2(i,1))/garums
      enddo
      read(*,*)
      deallocate(zm3)
#endif


      call getritzvectors(n_local,nsize,ndiv,Hx,blockH,extraHx,nstart,npw_local)
      call getritzvectors(n_local,nsize,ndiv,Sx,blockH,extraSx,nstart,npw_local)

      allocate(zvec(n_local))
      maxresid=0d0
      do i=1,ndiv
        zvec=extraHx(:,i)-rd(nstart-1+i)*extraSx(:,i)
        zsum=zdotc(npw,zvec(1),1,zvec(1),1)
#ifdef MPI
!        call MPI_ALLREDUCE(zsum, zsum2, 1,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!        zsum=zsum2
#endif
        if (nloall.gt.0) then
          residlen(i)=dble(zsum+zdotc(nloall,zvec(npw_local+1),1,zvec(npw_local+1),1))
        else
          residlen(i)=dble(zsum)
        endif
        maxresid=max(residlen(i),maxresid)
      enddo
      deallocate(zvec)
!     write(*,*) rd(nstart+1:nstart+ndiv)
      if (rank.eq.0) write(*,*) ndiv,maxresid,sum(rd(nstart:nstart+ndiv-1)),calls+1

      deallocate(BlockS,BlockH)
call timesec(time2)
#ifdef MEMORY_REPORT
if (rank.eq.0) then
 call memory_report
 write(*,*) 'Sx',sizeof(Sx)
 write(*,*) 'Hx',sizeof(Hx)
 write(*,*) 'trialvec',sizeof(trialvec)
 write(*,*) 'extraSx',sizeof(extraSx)
 write(*,*) 'extraHx',sizeof(extraHx)
 write(*,*) 'evecfv',sizeof(evecfv)
 write(*,*) 'apwi',sizeof(system%apwi)
endif
#endif

      calls=1
       newsum=sum(rd(nstart:nstart+ndiv-1))
       keepworking=.true.
       it=0

       do while (keepworking.and.(maxresid.gt.1d-20))
         it=it+1 


         ii=1
         do while (ii.lt.nblocks)
          ii=ii+1
          nadd=0
          do i=1,ndiv
             if (residlen(i).gt.1d-16) then
              nadd=nadd+1
              do j=1,npw
                trialvec(j,nadd+nsize)= (extraHx(j,i)-rd(nstart+i-1)*extraSx(j,i))/(hdiag(j)-rd(nstart+i-1)*sdiag(j))
              enddo
             endif
          enddo
      Hx(:,nsize+1:nsize+ndiv)=0d0
      Sx(:,nsize+1:nsize+ndiv)=0d0
 
!          nadd=ndiv
          evalfv(1:ndiv)=rd(nstart:nstart+ndiv-1)

if (nadd.ne.0) then
      nsize=nsize+nadd
      calls=calls+1
      allocate(BlockS(nsize,nsize))
      allocate(BlockH(nsize,nsize))
      call GSortho(n_local,nsize-nadd-nloall,nsize-nloall,trialvec(:,nloall+1:))
!      call GSortho(n,nsize-nadd,nsize,trialvec)

if (.true.) then
      call HapwSapw(n_local,npw,nadd,system,trialvec(:,nsize-nadd+1:nsize),Hx(:,nsize-nadd+1:nsize),Sx(:,nsize-nadd+1:nsize),.true.)
      if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
        call paceapw(n_local,npw,nadd,nstsv,0.25d0,pace(:,:,ik),trialvec(1:n_local,nsize-nadd+1:nsize) ,Hx(1:n_local,nsize-nadd+1:nsize))
      endif

      call innerproduct(n_local,nsize,nsize,trialvec,Hx,blockH,npw_local)
      call innerproduct(n_local,nsize,nsize,trialvec,Sx,blockS,npw_local)


      blockH(1:nsize-nadd,1:nsize-nadd)=extraH(1:nsize-nadd,1:nsize-nadd)
      do i=1,nsize-nadd
        do j=nsize-nadd+1,nsize
          blockH(j,i)=conjg(blockH(i,j))
        enddo
      enddo
      extraH(1:nsize,1:nsize)=BlockH(1:nsize,1:nsize)

      blockS(1:nsize-nadd,1:nsize-nadd)=extraS(1:nsize-nadd,1:nsize-nadd)
      do i=1,nsize-ndiv
        do j=nsize-ndiv+1,nsize
          blockS(j,i)=conjg(blockS(i,j))
        enddo
      enddo
      extraS(1:nsize,1:nsize)=BlockS(1:nsize,1:nsize)
else
!      call HapwSapw(n,npw,nsize,system,trialvec,Hx,Sx,.true.)
!      call innerproduct_lo(n,nsize,nsize,trialvec,Hx,blockH,npw,nloall)
!      call innerproduct_lo(n,nsize,nsize,trialvec,Sx,blockS,npw,nloall)


endif
      call diagH2(nsize,blockH,blockS,rd,info)
!if (info.eq.0) write(*,*) 'info=',info
!write(*,*) rd(1:ndiv)

if (info.eq.0) then
      call pickritzvectors(nsize,ndiv,rd,nstart)
      call getritzvectors(n_local,nsize,ndiv,trialvec,blockH,ritzvec,nstart,npw_local)
      call getritzvectors(n_local,nsize,ndiv,Hx,blockH,extraHx,nstart,npw_local)
      call getritzvectors(n_local,nsize,ndiv,Sx,blockH,extraSx,nstart,npw_local)


! length of residuals
      allocate(zvec(n_local))
      maxresid=0d0
      do i=1,ndiv
        zvec=extraHx(:,i)-rd(nstart-1+i)*extraSx(:,i)
        zsum=zdotc(npw,zvec(1),1,zvec(1),1)
#ifdef MPI
!        call MPI_ALLREDUCE(zsum, zsum2, 1,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!        zsum=zsum2
#endif
        if (nloall.gt.0) then
          residlen(i)=dble(zsum+zdotc(nloall,zvec(npw_local+1),1,zvec(npw_local+1),1))
        else
          residlen(i)=dble(zsum)
        endif
        maxresid=max(residlen(i),maxresid)
      enddo  
      deallocate(zvec) 



if (rank.eq.0)  write(*,*) nadd,maxresid,sum(rd(nstart:nstart-1+ndiv)),calls
!call barrier

#ifdef garums
      do i=1,nsize !nstart,ndiv+nstart
        garums=0d0
        summa=0d0
        do j=1,n-npw
          garums=garums+abs(blockH(j,i))**2
          summa=summa+abs(blockH(j,i))**2
        enddo
        do j=n-npw+1,n-npw+nusedsingular
          garums=garums+abs(blockH(j,i))**2
!          summa=summa+abs(blockH(j,i))**2
        enddo
        do j=n-npw+nusedsingular+1,nsize
          garums=garums+abs(blockH(j,i))**2
          summa=summa+abs(blockH(j,i))**2
        enddo
        write(*,*) i,summa/garums
       write(*,*) blockH(j,i)
      enddo
      read(*,*)
#endif

endif
      deallocate(BlockS,BlockH)


endif

          oldsum=newsum
          newsum=sum(rd(nstart:nstart+ndiv-1))
!          if (abs(newsum-oldsum).lt.input%groundstate%solver%epsarpack*nstfv) then
!          if ((nadd.eq.0).or.(calls.eq.80)) then
          if ((maxresid.lt.1d-16).or.(calls.eq.12)) then
            keepworking=.false.
            ii=nblocks
          endif
          if (info.ne.0) ii=nblocks
         enddo
    if (keepworking) then
       nsize=0
!       trialvec(:,nsize+1:nsize+nloall)=0d0
       trialvec=0d0
       nadd=0
       do i=1,n-npw

!         if (dble(hdiag(npw+i)).lt.5d0) then
           nadd=nadd+1
           trialvec(i+npw_local,nsize+nadd)=1d0
!         endif
       enddo
       nsize=nsize+nadd
!       call GSortho(n,0,nsize,trialvec)
       call HloSlo(n_local,npw_local,nsize,system,trialvec,Hx,Sx)      
       if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
        call pacelo(n_local,npw,nsize,nstsv,0.25d0,pace(:,:,ik),trialvec ,Hx )
       endif

       trialvec(1:npw_local,nsize+1:nsize+ndiv)=ritzvec(1:npw_local,1:ndiv)
       nsize=nsize+ndiv

      if (nusedsingular.ne.0) then
       trialvec(1:npw_local,nsize+1:nsize+nusedsingular)=singular(1:npw_local,1:nusedsingular,ik)
       nsize=nsize+nusedsingular
      endif

       call GSortho(n_local,0,nsize-nloall,trialvec(:,nloall+1:nsize))
       call HapwSapw(n_local,npw,nsize-nloall,system,trialvec(1:n_local,nloall+1:nsize),Hx(1:n_local,nloall+1:nsize),Sx(1:n_local,nloall+1:nsize),.true.)
       if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
         call paceapw(n_local,npw,nsize-nloall,nstsv,0.25d0,pace(:,:,ik),trialvec(1:n_local,nloall+1:nsize) ,Hx(1:n_local,nloall+1:nsize))
       endif

       allocate(BlockS(nsize,nsize))
       allocate(BlockH(nsize,nsize))
       extraH=0d0
       extraS=0d0
       call innerproduct(n_local,nsize,nsize,trialvec,Hx,blockH,npw_local)
       call innerproduct(n_local,nsize,nsize,trialvec,Sx,blockS,npw_local)
       extraH(1:nsize,1:nsize)=blockH(1:nsize,1:nsize)
       extraS(1:nsize,1:nsize)=blockS(1:nsize,1:nsize)
if (.false.) then
      call diagH2(nsize,blockH,blockS,rd,info)

      call pickritzvectors(nsize,ndiv,rd,nstart)
      call getritzvectors(n_local,nsize,ndiv,trialvec,blockH,ritzvec,nstart,npw_local)
      call getritzvectors(n_local,nsize,ndiv,Hx,blockH,extraHx,nstart,npw_local)
      call getritzvectors(n_local,nsize,ndiv,Sx,blockH,extraSx,nstart,npw_local)

      allocate(zvec(n_local))
      maxresid=0d0
      do i=1,ndiv
        zvec=extraHx(:,i)-rd(nstart-1+i)*extraSx(:,i)
        zsum=zdotc(npw,zvec(1),1,zvec(1),1)
#ifdef MPI
!        call MPI_ALLREDUCE(zsum, zsum2, 1,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!        zsum=zsum2
#endif
        residlen(i)=dble(zsum+zdotc(nloall,zvec(npw_local+1),1,zvec(npw_local+1),1))
        maxresid=max(residlen(i),maxresid)
      enddo
      deallocate(zvec)
      if (rank.eq.0) write(*,*) '*',ndiv,maxresid,sum(rd(1:ndiv))
endif
!      read(*,*)

       deallocate(BlockH,BlockS)
    endif     

       
  enddo
!stop
! Additional wavefunction filtering
if (.false.) then

      allocate(zm2(n,ndiv))
      allocate(idx(ndiv))
!      call Vx(n,npw,ndiv,system,ritzvec,zm2)
      j=0
      do i=1,ndiv
        write(*,*) nstart+i-1,rd(nstart+i-1),dble(zdotc(n,ritzvec(1,i),1,zm2(1,i),1))
        if (dble(zdotc(n,ritzvec(1,i),1,zm2(1,i),1)).lt.3d0*mine0) then
!          idx(i)=ndiv-(i-j)+1
          idx(ndiv-(i-j)+1)=i
          write(*,*) i, ndiv-(i-j)+1,dble(zdotc(n,ritzvec(1,i),1,zm2(1,i),1))
        else
          j=j+1
!          idx(i)=j
          idx(j)=i
        endif
!        write(*,*) nstart+i-1,rd(nstart+i-1),dble(zdotc(n,ritzvec(1,i),1,zm2(1,i),1))
      enddo
      if (j.eq.ndiv) then
        evalfv(1:ndiv)=rd(nstart:nstart+ndiv-1)
        evecfv (1:n, 1:nstfv)=ritzvec(1:n,1:nstfv)
      else
        write(*,*) 'new list of evals:'
        do i=1,j
          evalfv(i)=rd(nstart-1+idx(i))
          evecfv (1:n, i)=ritzvec(1:n,idx(i))
          write(*,*) evalfv(i), dble(zdotc(n,ritzvec(1,idx(i)),1,zm2(1,idx(i)),1))
        enddo
        do i=j+1,ndiv
          evalfv(i)=1d5
          evecfv (1:n, i)=ritzvec(1:n,idx(i))
          write(*,*) evalfv(i)
        enddo

      endif
      deallocate(idx,zm2)
else
      evalfv(1:ndiv)=rd(nstart:nstart+ndiv-1)
!write(*,*)
!do i=1,ndiv
!  write(*,*) evalfv(i)
!enddo
!      if (procs_per_kpt.eq.1) then
        evecfv (1:npw, 1:nstfv)=ritzvec(1:npw,1:nstfv)
        if (nloall.ne.0) then
          evecfv (npw+1:n, 1:nstfv)=ritzvec(n_local-nloall+1:n_local,1:nstfv)
        endif
!      else
#ifdef MPI
#ifdef TIMINGS
!        call timesec(time1)
#endif
!        allocate(bufrecv(npw_local*procs_per_kpt))
!        do i=1,nstfv
!          call MPI_GATHER(ritzvec(1,i), npw_local, MPI_DOUBLE_COMPLEX, bufrecv, npw_local, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr)
!          if (rank.eq.0) evecfv(1:npw,i)=bufrecv(1:npw)
!        enddo
!        deallocate(bufrecv)
#ifdef TIMINGS
!        call timesec(time2)
!        if (rank_in_kpt.eq.0) write(*,*) 'collecting WFs with MPI', time2-time1
#endif
!        if (rank_in_kpt.eq.0) evecfv(npw+1:n,1:nstfv)=ritzvec(npw_local+1:n_local,1:nstfv)
#endif
!      endif
!write(*,*) 'nsize',nsize
endif

      write(*,*) '***** itnum=',it,calls
      call timesec(tsb)
      write(*,*) 'iterations',tsb-tsa
      timefv=timefv+tsb-tsa
      Return
End Subroutine davidson 

