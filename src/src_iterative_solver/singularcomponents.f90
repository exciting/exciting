Subroutine singularcomponents (system, nst, evecfv, evalfv,ik)
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
      use modxs, only : fftmap_type

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
      integer, intent(in) :: ik
!
  ! local variables
      Logical :: packed
      Type (HermitianMatrix) :: fm
!
      Integer :: n
      Real(8) :: cpu0, cpu1, cpu2, tsa,tsb
      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
  !IO vars
      Integer :: koffset, recl
      Character (256) :: outfilenamestring, filetag
      External outfilenamestring
!
  !ARPACK Interface vars
      Integer :: ido, nev, ncv, lworkl, info, infoznaupd, info2, j, i,it
      Integer :: nevmax, ncvmax, nmax
      Integer :: nconv, maxitr, ishfts, mode, ldv
      Integer :: iparam (11), ipntr (14)
      Complex (8), Allocatable :: resid (:), v (:, :), workev (:), &
     & workl (:), d (:), buf(:)
      Complex (8), Pointer :: workd (:)
      Real (8), Allocatable :: rwork (:), rd (:)
      Integer, Allocatable :: idx (:)
      Complex (8) :: sigma
      Character :: bmat * 1, which * 2
      Real (8) :: tol
      Logical :: rvec
      Logical :: select (nmat(1, ik))
      Complex (8), Pointer :: vin (:), vout (:)

      logical :: ImproveInverse
      type (HermitianMatrix) :: zm,cc
      Complex (8), Allocatable :: bres(:),vecupd(:)
      integer :: ii,ib,nounter
      real(8) :: berr
      
      real(8) ta,tb

      real(8) :: DegeneracyThr
      integer :: blockstart,blocksize
      Complex (8), Allocatable :: h(:,:),blockH(:,:),blockS(:,:),Hx(:,:),Sx(:,:)
      Complex (8), Allocatable :: zm2(:,:),diag(:),zm3(:,:),zm4(:,:)
      Complex (8), allocatable :: trialvec(:,:),sdiag(:),resvec(:)
      real(8), allocatable :: evals(:),sortedevals(:),residlen(:)
      logical :: SeparateDegenerates,sorted

      real(8), allocatable :: rdiag(:)
      complex(8), allocatable :: cdiag(:),cinvdiag(:),zvec(:),rrvec(:),pvec(:),qvec(:),extraS(:,:),extraH(:,:),extraSx(:,:),extraHx(:,:),ritzvec(:,:)
      complex(8), allocatable :: zfftcf(:)
      Complex (8), Allocatable :: work (:),cwork(:), zax(:,:) !,singular(:,:)
      integer, allocatable :: iwork(:)

      integer :: counter
      character(len=4) :: strInvertMethod

      logical :: gev,keepworking
      integer :: method

      Type (HermitianMatrix), pointer ::pinverse
! sandbox variables :)
      integer :: ig
      real(8) :: summa,maxdiff
      complex(8) :: rho,beta,alpha,rhoold,zsum,zsum2
      logical :: usecg, calcsingular
 
      integer :: lwork
      complex(8) :: worksize

      real(8) :: oldsum,newsum,norm,sft,newrd,maxkin,potl,time1,time2,res
      integer :: ndiv,nblocks,calls,npw,is,ia,blo,bhi,l,lm,ias,nsize,nadd,n_local

      complex(8), external :: zdotc
      integer, allocatable :: offset(:),ngklist(:),ibuf(:)

      type(fftmap_type) :: fftmap



call timesec(tsa)

      tol = input%groundstate%solver%epsarpack
      packed=.false.

      npw=ngk(1,ik)

      n=nmat(1,ik)
      n_local=npw


if (nsingular.eq.-1) then
if (.true.) then
  nsingular=0
  do is=1,nspecies
    nsingular=nsingular+natoms(is)*int(1.2d0*2.8d-3*(rmt(is)*gkmax)**4)
  enddo
else
  read(*,*) nsingular
endif
if (nsingular.eq.0) nsingular=1
endif
ndiv=nsingular

      nblocks=8 !10

if (.not.(allocated(singular))) then
 calcsingular=.true.
else
 calcsingular=(sum(singular(:,:,ik)).eq.0d0)
endif
if  (calcsingular) then
!if (.true.) then

allocate(sdiag(n_local))

sdiag= 0d0
do i=1,npw
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      sdiag(i)=sdiag(i)+zdotc(mt_hscf%maxaa,system%apwi(1,i,ias),1,system%apwi(1,i,ias),1)
    enddo
  enddo
enddo

      sdiag(1:npw)=sdiag(1:npw)+ cfunig0
      allocate(residlen(nblocks*ndiv))
      allocate(rd(nblocks*ndiv))
      allocate(evals(ndiv))
      allocate(trialvec(n_local,nblocks*ndiv))
      allocate(extraS(nblocks*ndiv,nblocks*ndiv))
      allocate(extraH(nblocks*ndiv,nblocks*ndiv))
      allocate(ritzvec(n_local,nblocks*ndiv))
      allocate(resvec(n_local))


      if (associated(system%overlap%za)) then 
        allocate(zfftcf(1))
        nullify(fftmap%igfft)
      else
        call genfftmap(fftmap,2.01d0*gkmax) 
        allocate(zfftcf(fftmap%ngrtot))
        zfftcf=0d0
        do ig=1, fftmap%ngvec
          zfftcf(fftmap%igfft(ig))=cfunig(ig)
        end do
        call zfftifc(3, fftmap%ngrid, 1, zfftcf)
      endif
      

      extraH=zero
      extraS=zero
      ritzvec=zero
      trialvec=zero

      

      do i=1,ndiv
        trialvec(i,i)=one
        trialvec(i+1,i)=0.5d0
        trialvec(i+1,i)=0.25d0
      enddo

      allocate(zm2(n_local,nblocks*ndiv))
      allocate(Sx(n_local,nblocks*ndiv))
      Sx=0d0
      calls=1
      nsize=ndiv

      allocate(BlockS(nsize,nsize))
      allocate(BlockH(nsize,nsize))

      call GSortho(n_local,0,ndiv,trialvec) 

      call OverlapX(n_local,npw,ndiv,system,fftmap,zfftcf,trialvec,Sx,.true.)
!write(*,*) 'sumsx',sum(Sx)
!stop
      call innerproduct(n_local,nsize,nsize,trialvec,Sx,blockH)
      call innerproduct(n_local,nsize,nsize,trialvec,trialvec,blockS)
      extraH(1:ndiv,1:ndiv)=blockH(1:ndiv,1:ndiv)
      extraS(1:ndiv,1:ndiv)=blockS(1:ndiv,1:ndiv)
!if (rank.eq.0) then
      call diagH(ndiv,blockH,blockS,rd,info)
!endif
#ifdef MPI
!      call MPI_BCAST(blockH, 2*ndiv**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
!      call MPI_BCAST(rd, ndiv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
#endif

      call getritzvectors(n_local,nsize,nsize,trialvec,blockH,ritzvec,1,n_local)
      call getritzvectors(n_local,nsize,nsize,Sx,blockH,zm2,1,n_local)
      deallocate(BlockS,BlockH)

       evals(1:ndiv)=rd(1:ndiv)
       newsum=sum(rd(1:ndiv))
write(*,*) newsum
#ifdef MEMORY_REPORT
if (rank.eq.0) then
 call memory_report
 write(*,*) 'Sx',sizeof(Sx)
 write(*,*) 'zm2',sizeof(zm2)
 write(*,*) 'trialvec',sizeof(trialvec)
 write(*,*) 'apwi',sizeof(system%apwi)
endif
#endif

       keepworking=.true.
       it=0

       do while (keepworking)

         it=it+1
         ii=1
         do while (ii.lt.nblocks)
          ii=ii+1
          do i=1,ndiv
            do j=1,npw
             trialvec(j,i+(ii-1)*ndiv)=(zm2(j,i)-rd(i)*ritzvec(j,i))/(sdiag(j)-rd(i))
            enddo
          enddo
          calls=calls+1
          nsize=ii*ndiv
          allocate(BlockS(nsize,nsize))
          allocate(BlockH(nsize,nsize))
          allocate(zm3(n_local,nsize))
!!write(*,*) it,ii,'res',sum(trialvec)
          call GSortho(n_local,nsize-ndiv,nsize,trialvec) 
!!write(*,*) it,ii,'GSortho',sum(trialvec)
          call OverlapX(n_local,npw,ndiv,system,fftmap,zfftcf,trialvec(:,nsize-ndiv+1:nsize),zm3,.true.)
!!write(*,*) it,ii,'sumsx',sum(Sx)
blockH=0d0
          call innerproduct(n_local,nsize,ndiv,trialvec,zm3,blockH(:,nsize-ndiv+1:nsize))

          blockH(1:nsize-ndiv,1:nsize-ndiv)=extraH(1:nsize-ndiv,1:nsize-ndiv)

          do i=1,nsize-ndiv
            do j=nsize-ndiv+1,nsize
              blockH(j,i)=conjg(blockH(i,j))
            enddo
          enddo
          Sx(1:n_local,nsize-ndiv+1:nsize)=zm3(1:n_local,1:ndiv)
          extraH(1:nsize,1:nsize)=BlockH(1:nsize,1:nsize)

          call innerproduct(n_local,nsize,ndiv,trialvec,trialvec(:,nsize-ndiv+1:nsize),blockS(:,nsize-ndiv+1:nsize))
          blockS(1:nsize-ndiv,1:nsize-ndiv)=extraS(1:nsize-ndiv,1:nsize-ndiv)
          do i=1,nsize-ndiv
            do j=nsize-ndiv+1,nsize
              blockS(j,i)=conjg(blockS(i,j))
            enddo
          enddo
          extraS(1:nsize,1:nsize)=BlockS(1:nsize,1:nsize)

!          if (rank.eq.0) then
            call diagH(nsize,blockH,blockS,rd,info)
!          endif
!#ifdef MPI
!          call MPI_BCAST(blockH, 2*nsize**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(rd, nsize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(info, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
!#endif
!!write(*,*) it,ii,'diagh',sum(rd)
          if (info.eq.0) then
            call getritzvectors(n_local,nsize,nsize,trialvec,blockH,ritzvec,1,n_local)
            call getritzvectors(n_local,nsize,nsize,Sx,blockH,zm2,1,n_local)
!trialvec(gkhi_local+1:nlocal,
!           call getritzvectors2(n_local,nsize,ndiv,trialvec,blockH,ritzvec,nstart,npw_local)
!            Sx(1:n_local,1:nsize)=zm2(1:n_local,1:nsize)
          endif
          deallocate(zm3,BlockS,BlockH)

!          if ((abs((sum(evals(1:ndiv))-sum(rd(1:ndiv)))/sum(rd(1:ndiv))).lt.1d-2).and. &
!              (abs((evals(1)-rd(1))/rd(1)).lt.1d-2)) then
!            ii=nblocks
!            keepworking=.false.
!          endif

          res=0d0
          do i=1,ndiv
            resvec(1:npw)=(zm2(1:npw,i)-rd(i)*ritzvec(1:npw,i))
            zsum=zdotc(npw,resvec,1,resvec,1)
!#ifdef MPI
!            call MPI_ALLREDUCE(zsum, zsum2, 1,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!            zsum=zsum2
!#endif
            residlen(i)=dble(zsum)
            if (dble(zsum).gt.res) res=dble(zsum)
!            write(*,*) 'res',i,rd(i),residlen(i)
          enddo
!read(*,*)

          keepworking=(abs(res).gt.1d-8)
          write(*,*) sum(rd(1:ndiv)),res,nsize,ii

          evals(1:ndiv)=rd(1:ndiv)
          if ((info.ne.0).or.(.not.keepworking)) ii=nblocks
          
         enddo


!         nsize=ndiv
         trialvec=zero
         trialvec(1:n_local,1:ndiv)=ritzvec(1:n_local,1:ndiv)
         call OverlapX(n_local,npw,ndiv,system,fftmap,zfftcf,trialvec,zm2,.true.)
         Sx(1:n_local,1:ndiv)=zm2(1:n_local,1:ndiv)
         allocate(BlockH(ndiv,ndiv))
         call innerproduct(n_local,ndiv,ndiv,trialvec,zm2,blockH(1:ndiv,1:ndiv))

         extraH=zero
         extraH(1:ndiv,1:ndiv)=blockH(1:ndiv,1:ndiv)
         extraS=zero
         do i=1,ndiv
           extraS(i,i)=one
         enddo
         deallocate(blockH)

       enddo
!if (ik.ne.1) then
! write(*,*) 'fix the size of singular(...)'
! stop
!endif
write(*,*) calls,' calls'
if (.not.(allocated(singular))) then
  i=ndiv
  do while ((i.gt.1).and.(rd(i).gt.3d-1))
    i=i-1
  enddo
!  i=1
  if (nsingular.ne.0) then
    nsingular=i
  else
    nsingular=1
  endif
  allocate(singular(ngkmax,nsingular,nkpt))
  allocate(evalsingular(nsingular,nkpt))

  singular=zzero
endif
!stop

deallocate(zm2)
if (associated(fftmap%igfft)) deallocate(fftmap%igfft)
deallocate(zfftcf)

singular(1:n_local,1:nsingular,ik)=ritzvec(1:n_local,1:nsingular)
evalsingular(1:nsingular,ik)=rd(1:nsingular)

write(*,*) nsingular, 'singular components'
write(*,*) '***********'
write(*,*) 'highest eigenvalue among singular components',rd(nsingular)
write(*,*) rd(1:nsingular)
deallocate(sdiag)
!read(*,*)
endif

      call timesec(tsb)
      write(*,*) 'singular components',tsb-tsa
      timefv=timefv+tsb-tsa
      Return
End Subroutine singularcomponents 

