! Copyright (C) 2015-2023 exciting team (Berlin and Riga)
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine davidson (system, nst, evecfv, evalfv,ik)
!> The subroutine calculates the required number of the smallest eigenpairs.
!> The calculation is performed using a modified Davidson algorithm which requires that the operations H|psi> and S|psi> are defined.  
!> If the data structure system contains the Hamiltonian and overlap matrices, H|psi> and S|psi> are performed explicitly as matrix-vector multiplication.
!> If the matrices are not given, they are never constructed explicitly. 
!> The algorithm constructs the subspace that contains the smallest eigenvectors of S, all local orbitals and an expanding set of (L)APW degrees of freedom. 
!> 
      Use constants, Only : zzero,zone 
      Use modinput
      Use modfvsystem
      Use modmpi, only: terminate_mpi_env, mpiglobal
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
      Use mod_kpoint, Only: nkpt
      use modxs, only : fftmap_type
      Implicit None
! arguments
      Real (8), Intent (out) :: evalfv (nst)                    ! eigenvalues
      Complex (8), Intent (inout) :: evecfv (nmatmax, nst)      ! eigenvectors
      Type (evsystem) :: system                                 ! a datastructure with either the Hamiltonian and overlap matrices or matching coefficients. 
      integer, intent(in) :: nst                                ! number of eigenpairs required
      integer, intent(in) :: ik                                 ! which k-vector is considered
!
! local variables
!
      Integer :: n
      Real(8) :: tsa,tsb
      integer :: j,i,nstart,info
      Real (8), Allocatable :: rd (:),residlen(:)
      Real (8), parameter :: tol=1d-16
      Real (8) :: maxresid

      integer :: ii
      
      Complex (8), Allocatable :: blockH(:,:),blockS(:,:),Hx(:,:),Sx(:,:)
      Complex (8), allocatable :: trialvec(:,:)

      complex(8), allocatable :: zvec(:),extraS(:,:),extraH(:,:),extraSx(:,:),extraHx(:,:),ritzvec(:,:)
      Complex (8), Allocatable :: sdiag(:),hdiag(:) 
      integer :: ig
      complex(8) :: zsum
 
      real(8) :: oldsum,newsum,time1,time2
      integer :: ndiv,nblocks,calls,npw,is,ia,ias,nsize,nadd,ilo,m,if1,nloall,n_local,npw_local,nusedsingular,l
      complex(8), external :: zdotc
      complex(8), allocatable :: zfftcf(:),zfftveff(:),zfftmeff(:)
      type(fftmap_type) :: fftmap


      if ((input%groundstate%outputlevel.eq."high").and.(mpiglobal%rank.eq.0))  write(*,*) 'ik=',ik
      call timesec(tsa)

! Initialise block sizes
      npw=ngk(1,ik)
      n=nmat(1,ik)   
      nloall=n-npw
      npw_local=npw 
      n_local=nloall+npw_local

! Diagonal elements of the Hamiltonian and overlap matrices
      allocate(sdiag(n_local))
      allocate(hdiag(n_local))
      sdiag=0d0
      hdiag=0d0

! Was the matrix constructed
      if (associated(system%hamilton%za)) then 
! Yes, then it is trivial
        do i=1,n
          hdiag(i)=system%hamilton%za(i,i)
          sdiag(i)=system%overlap%za(i,i)
        enddo
      else
! No, let's assemble the diagonal then
! First, the overlap
        sdiag(1:npw)= cfunig(1)

        Do i=1,npw
          Do is = 1, nspecies
            Do ia = 1, natoms (is)
              ias = idxas (ia, is)
              sdiag(i)=sdiag(i)+zdotc(mt_hscf%maxaa,system%apwi(1,i,ias),1,system%apwi(1,i,ias),1)
            enddo
          enddo
        enddo

        sdiag(npw_local+1:n_local)=1d0

! Now the Hamiltionian
        hdiag(1:npw)= veffig(1)

        if (input%groundstate%ValenceRelativity.ne."none") then
          do i=1,npw
            hdiag(i)=hdiag(i)+0.5d0*dot_product (current_vgkc(:, i), current_vgkc(:, i))*meffig(1)
          enddo
        else
          do i=1,npw
            hdiag(i)=hdiag(i)+0.5d0*dot_product (current_vgkc(:, i), current_vgkc(:, i))*cfunig(1)
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
              call zgemv('N',mt_hscf%maxaa,mt_hscf%maxaa,zone,mt_hscf%main%aa(1,1,ias),mt_hscf%maxaa,system%apwi(1,i,ias),1,zzero,zvec,1)
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

! How will we apply S|psi> and H|psi>?
      if (associated(system%overlap%za)) then
! We have the matrices, no need to initialise FFT
        allocate(zfftcf(1))
        allocate(zfftveff(1))
        allocate(zfftmeff(1))
        nullify(fftmap%igfft)
      else
! Initialise a new FFT grid assuming gmax = 2.01 gkmax.
        call genfftmap(fftmap,2.01d0*gkmax)

        allocate(zfftcf(fftmap%ngrtot))
        zfftcf=0d0
        do ig=1, fftmap%ngvec
          zfftcf(fftmap%igfft(ig))=cfunig(ig)
        end do
        call zfftifc(3, fftmap%ngrid, 1, zfftcf)

        allocate(zfftveff(fftmap%ngrtot))
        zfftveff=0d0
        do ig=1, fftmap%ngvec
          zfftveff(fftmap%igfft(ig))=veffig(ig)
        end do
        call zfftifc(3, fftmap%ngrid, 1, zfftveff)

        if (input%groundstate%ValenceRelativity.ne."none") then
          allocate(zfftmeff(fftmap%ngrtot))
          zfftmeff=0d0
          do ig=1, fftmap%ngvec
            zfftmeff(fftmap%igfft(ig))=meffig(ig)
          end do
          call zfftifc(3, fftmap%ngrid, 1, zfftmeff)
        else
          allocate(zfftmeff(1))
        endif

      endif




      nusedsingular=nsingular 
      ndiv=nst
      nblocks=12 ! Expand the subspace up to 12 times
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

! Initialise the subspace 
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
! The following bit will be necessary for supporting ACE
!       if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
!        call pacelo(n_local,npw_local,nsize,nstsv,input%groundstate%hybrid%excoeff,pace(:,:,ik),trialvec ,Hx )
!       endif

! -> 2. singular components
       if (nusedsingular.ne.0) then
         trialvec(1:npw_local,nsize+1:nsize+nusedsingular)=singular(1:npw_local,1:nusedsingular,ik)
         nsize=nsize+nusedsingular
         call GSortho(n_local,0,nsize-nloall,trialvec(:,nloall+1:))
       endif
! -> 3. initial guess for wavefunctions
       if (dble(sum(evecfv)).ne.0d0) then ! check whether the old eigenvectors exist
!       wavefunctions from previous scf iterations
         trialvec(1:npw,nsize+1:nsize+nst)=evecfv(1:npw, 1:nst)
         nsize=nsize+nst
       else
!       or just a guess
         if (4d0*nst.ge.npw) then
           write(*,*) 'The number of required eigenvalues is too high for the davidson eigensolver'
           write(*,*) 'The allowed maximum is one fourth of the number of LAPWs.'
           call terminate_mpi_env(mpiglobal)
         endif
         do i=1,ndiv
           trialvec((i-1)*4+1,nsize+i)=zone
           trialvec((i-1)*4+1+1,nsize+i)=zone !0.5d0
           trialvec((i-1)*4+1+2,nsize+i)=zone !0.25d0
           trialvec((i-1)*4+1+3,nsize+i)=zone !0.125d0
         enddo
         nsize=nsize+ndiv
       endif

! Construct the (Ritz) eigenproblem within the subspace 

call timesec(time1)

      allocate(BlockS(nsize,nsize))
      allocate(BlockH(nsize,nsize))
      BlockH=zzero
      BlockS=zzero

      call GSortho(n_local,nsize-ndiv-nloall,nsize-nloall,trialvec(:,nloall+1:))
      call HapwSapw(n_local,npw,nsize-nloall,system,fftmap,zfftcf,zfftveff,zfftmeff,trialvec(:,nloall+1:nsize),Hx(:,nloall+1:nsize),Sx(:,nloall+1:nsize),.true.)

! The following bit will be necessary for supporting ACE
!      if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
!        call paceapw(n_local,npw,nsize-nloall,nstsv,input%groundstate%hybrid%excoeff,pace(:,:,ik),trialvec(1:n_local,nloall+1:nsize) ,Hx(1:n_local,nloall+1:nsize))
!      endif

      call innerproduct(n_local,nsize,nsize,trialvec,Hx(:,1:nsize),blockH(:,1:nsize))
      call innerproduct(n_local,nsize,nsize,trialvec,Sx(:,1:nsize),blockS(:,1:nsize))
      extraH(1:nsize,1:nsize)=BlockH(1:nsize,1:nsize)
      extraS(1:nsize,1:nsize)=BlockS(1:nsize,1:nsize)

! The initial subspace is set up
! Now diagonalise and get the current estimate of the eigenvectors 

      call diagH2(nsize,blockH,blockS,rd,info)
      call pickritzvectors(nsize,ndiv,rd,nstart)
      call getritzvectors(n_local,nsize,ndiv,trialvec,blockH,ritzvec,nstart,npw_local)

      call getritzvectors(n_local,nsize,ndiv,Hx,blockH,extraHx,nstart,npw_local)
      call getritzvectors(n_local,nsize,ndiv,Sx,blockH,extraSx,nstart,npw_local)

! Calculate the residuals

      allocate(zvec(n_local))
      maxresid=0d0
      do i=1,ndiv
        zvec=extraHx(:,i)-rd(nstart-1+i)*extraSx(:,i)
        zsum=zdotc(npw,zvec(1),1,zvec(1),1)
        if (nloall.gt.0) then
          residlen(i)=dble(zsum+zdotc(nloall,zvec(npw_local+1),1,zvec(npw_local+1),1))
        else
          residlen(i)=dble(zsum)
        endif
        maxresid=max(residlen(i),maxresid)
      enddo
      deallocate(zvec)
      deallocate(BlockS,BlockH)

      if ((input%groundstate%outputlevel.eq."high").and.(mpiglobal%rank.eq.0)) write(*,*) ndiv,maxresid,sum(rd(nstart:nstart+ndiv-1)),calls+1

call timesec(time2)

      calls=1
      newsum=sum(rd(nstart:nstart+ndiv-1))

! The Davidson iterations
! Expand the subspace up to nblocks times according to the Davidson's recipe assuming diagonally dominant Hamiltonian
      ii=0
      do while (ii.lt.nblocks)
        ii=ii+1
        nadd=0
        do i=1,ndiv
          if (residlen(i).gt.tol) then ! If the residual is small enough we skip the the eigenvector
            nadd=nadd+1
            do j=1,npw
              trialvec(j,nadd+nsize)= (extraHx(j,i)-rd(nstart+i-1)*extraSx(j,i))/(hdiag(j)-rd(nstart+i-1)*sdiag(j))
            enddo
          endif
        enddo

        evalfv(1:ndiv)=rd(nstart:nstart+ndiv-1)


        Hx(:,nsize+1:nsize+ndiv)=0d0
        Sx(:,nsize+1:nsize+ndiv)=0d0

        if (nadd.ne.0) then
          nsize=nsize+nadd
          calls=calls+1
          allocate(BlockS(nsize,nsize))
          allocate(BlockH(nsize,nsize))
          call GSortho(n_local,nsize-nadd-nloall,nsize-nloall,trialvec(:,nloall+1:))
          call HapwSapw(n_local,npw,nadd,system,fftmap,zfftcf,zfftveff,zfftmeff,trialvec(:,nsize-nadd+1:nsize),Hx(:,nsize-nadd+1:nsize),Sx(:,nsize-nadd+1:nsize),.true.)
! The following bit will be necessary for supporting ACE 
!      if (allocated(pace).and.(.not.associated(system%hamilton%za))) then
!        call paceapw(n_local,npw,nadd,nstsv,input%groundstate%hybrid%excoeff,pace(:,:,ik),trialvec(1:n_local,nsize-nadd+1:nsize) ,Hx(1:n_local,nsize-nadd+1:nsize))
!      endif

          call innerproduct(n_local,nsize,nsize,trialvec,Hx,blockH)
          call innerproduct(n_local,nsize,nsize,trialvec,Sx,blockS)


          blockH(1:nsize-nadd,1:nsize-nadd)=extraH(1:nsize-nadd,1:nsize-nadd)
          extraH(1:nsize,1:nsize)=BlockH(1:nsize,1:nsize)

          blockS(1:nsize-nadd,1:nsize-nadd)=extraS(1:nsize-nadd,1:nsize-nadd)
          extraS(1:nsize,1:nsize)=BlockS(1:nsize,1:nsize)

          call diagH2(nsize,blockH,blockS,rd,info)

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
              if (nloall.gt.0) then
                residlen(i)=dble(zsum+zdotc(nloall,zvec(npw_local+1),1,zvec(npw_local+1),1))
              else
                residlen(i)=dble(zsum)
              endif
              maxresid=max(residlen(i),maxresid)
            enddo  
            deallocate(zvec) 

            if ((input%groundstate%outputlevel.eq."high").and.(mpiglobal%rank.eq.0))  write(*,*) nadd,maxresid,sum(rd(nstart:nstart-1+ndiv)),calls

          elseif (mpiglobal%rank.eq.0) then
            write(*,*) 'Subspace diagonalisation failed in davidson.f90.'
            write(*,*) 'Is the subspace linearly dependent?'
            write(*,*) 'info=',info
            call terminate_mpi_env(mpiglobal) 
          endif
          deallocate(BlockS,BlockH)
        endif

        oldsum=newsum
        newsum=sum(rd(nstart:nstart+ndiv-1))

        if ((maxresid.lt.tol).or.(calls.eq.nblocks)) then
          ii=nblocks
        endif
       enddo

! Cleaning up
      if (associated(fftmap%igfft)) deallocate(fftmap%igfft)
      deallocate(zfftcf)
      deallocate(zfftveff)

      evalfv(1:ndiv)=rd(nstart:nstart+ndiv-1)
      evecfv (1:npw, 1:nstfv)=ritzvec(1:npw,1:nstfv)
      if (nloall.ne.0) then
        evecfv (npw+1:n, 1:nstfv)=ritzvec(n_local-nloall+1:n_local,1:nstfv)
      endif

      call timesec(tsb)

      timefv=timefv+tsb-tsa
      if ((input%groundstate%outputlevel.eq."high").and.(mpiglobal%rank.eq.0)) then
        write(*,*) '***** HapwSapw calls=',calls
        write(*,*) 'iterations',tsb-tsa
      endif
      Return
End Subroutine davidson 


