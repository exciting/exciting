!> The subroutine calculates (i) the smallest eigenvalues of the APW-APW block of the overlap matrix 
!> and (ii) the corresponding eigenvectors. 
!> The calculation is performed using the Davidson algorithm which requires that the operation S|x> is defined.  
!> If the data structure system contains the overlap matrix, S|x> is performed explicitly as matrix-vector multiplication.
!> If the overlap matrix is missing, it is never constructed explicitly. 
!> The outcome of the calculation is the eigenvectors that correspond to eigenvalues under 0.3 (hard-coded threshold). 
!> These eigenvectors are then stored in the global variable singular.
subroutine singularcomponents ( mpi_env, system, ik )
  use precision, only: dp, i32
  use constants, only : zzero, zone
  use modinput, only: input
  use modfvsystem, only: evsystem
  use exciting_mpi, only: mpiinfo
  use mod_atoms, only: natoms, nspecies, idxas
  use mod_eigensystem, only: nmat, singular, nsingular, mt_hscf
  use mod_timing, only: timefv
  use mod_Gkvector, only: gkmax, ngk, ngkmax
  use mod_Gvector, only: cfunig
  use mod_muffin_tin, only: rmt
  use mod_atoms, only: natmtot
  use mod_spin, only: nspnfv
  use mod_APW_LO, only: apwordmax
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_kpoint, only: nkpt
  use modxs, only : fftmap_type
  use m_zfftifc, only: zfftifc

  implicit none

  type(mpiinfo), intent(in) :: mpi_env
  type (evsystem) :: system
  integer, intent(in) :: ik

  real(dp) :: tsa, tsb, newsum, res
  integer(i32) :: ii, page, i, j, it, info, n, ig, ndiv, nblocks, calls, npw, is, ia, ias, nsize, n_local
  complex(dp), allocatable :: trialvec(:,:), sdiag(:), resvec(:), blockH(:,:), &
    blockS(:,:), Sx(:,:), zm2(:,:), zm3(:,:), extraS(:,:), extraH(:,:), ritzvec(:,:), &
    zfftcf(:)
  real(dp), allocatable :: evals(:), rd (:)
  logical :: keepworking, calcsingular
  complex(dp) :: zsum
  complex(dp), external :: zdotc
  type(fftmap_type) :: fftmap

  call timesec(tsa)

  npw=ngk(1,ik)
  n=nmat(1,ik)
  n_local=npw
  page = ik
  if ( allocated( singular ) ) then
    if ( size( singular, 3 ) == 1 ) page = 1
  endif

  ! How many singular components do we expect to get?
  ! The formula is just an empirical estimate. 
  if (nsingular.eq.-1) then
    nsingular=0
    do is=1,nspecies
      nsingular=nsingular+natoms(is)*int(1.2d0*2.8d-3*(rmt(is)*gkmax)**4)
    enddo
    if (nsingular.eq.0) nsingular=1
  endif
  ndiv=nsingular

  ! The subspace will be expanded nblocks times before the restart
  nblocks=8 

  if (.not.(allocated(singular))) then
    calcsingular=.true.
  else
    calcsingular=(sum(singular(:,:,page)).eq.0d0)
  endif

  ! skip the whole thing if the singular component were already calculated
  if (calcsingular) then

    allocate(sdiag(n_local))

    ! We need to store the diagonal of the overlap matrix (only the APW-APW part).
    ! We initialise it below regardless whether the overlap is constructed or not.
    sdiag= 0d0
    do i=1,npw
      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          sdiag(i)=sdiag(i)+zdotc(mt_hscf%maxaa,system%apwi(1,i,ias),1,system%apwi(1,i,ias),1)
        enddo
      enddo
    enddo

    sdiag(1:npw)=sdiag(1:npw)+ cfunig(1)

    allocate(rd(nblocks*ndiv))
    allocate(evals(ndiv))
    allocate(resvec(n_local))

    ! If the overlap matrix is not available, 
    ! initialise FFT grids with gmax = 2*gkmax.
    if (associated(system%overlap%za)) then 
      allocate(zfftcf(1))
      nullify(fftmap%igfft)
    else
      call genfftmap(fftmap,2.01_dp*gkmax) 
      allocate(zfftcf(fftmap%ngrtot))
      zfftcf=0d0
      do ig=1, fftmap%ngvec
        zfftcf(fftmap%igfft(ig))=cfunig(ig)
      end do
      call zfftifc(3, fftmap%ngrid, 1, zfftcf)
    endif
      
    ! Initialise arrays:
    ! trialvec - subspace
    ! extraH = <a|S|b>
    ! extraS = <a|b>
    ! ritzvec - eigenvectors of the Ritz problem extraH*x = sigma*extraS*x
    allocate(trialvec(n_local,nblocks*ndiv), source = zzero)
    allocate(extraS(nblocks*ndiv,nblocks*ndiv), source = zzero)
    allocate(extraH(nblocks*ndiv,nblocks*ndiv), source = zzero)
    allocate(ritzvec(n_local,nblocks*ndiv), source = zzero)
      
    ! Initial guess
    do i=1,ndiv
      trialvec(i,i)=zone
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

    ! calculate S|trialvec>
    call OverlapX(n_local,npw,ndiv,system,fftmap,zfftcf,trialvec,Sx,.true.)

    ! calculate <trialvec'|S|trialvec> and <trialvec'|trialvec>
    call innerproduct(n_local,nsize,nsize,trialvec,Sx,blockH)
    call innerproduct(n_local,nsize,nsize,trialvec,trialvec,blockS)

    ! store blockH and blockS, because diagonalisation will ruin it
    extraH(1:ndiv,1:ndiv)=blockH(1:ndiv,1:ndiv)
    extraS(1:ndiv,1:ndiv)=blockS(1:ndiv,1:ndiv)

    ! subspace diagonalisation blockH|x> = sigma*blockS|x>
    call diagH2(ndiv,blockH,blockS,rd,info)

    ! update |trialvec> 
    call getritzvectors(n_local,nsize,nsize,trialvec,blockH,ritzvec,1,n_local)
    call getritzvectors(n_local,nsize,nsize,Sx,blockH,zm2,1,n_local)
    deallocate(BlockS,BlockH)

    evals(1:ndiv)=rd(1:ndiv)
    newsum=sum(rd(1:ndiv))
    keepworking=.true.
    it=0

    ! The Davidson loop with restarts
    do while (keepworking)
      it=it+1
      ii=1
      ! The inner loop - expand the subspace and diagonalise
      do while (ii.lt.nblocks)
        ii=ii+1
        ! Expand the subspace
        ! We expect a diagonally dominant overlap, hence the formula
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
        ! Orthogonalise the new part of |trialvec> 
        ! Calculate S|trialvec> - yet again for the new part only
        call GSortho(n_local,nsize-ndiv,nsize,trialvec) 
        call OverlapX(n_local,npw,ndiv,system,fftmap,zfftcf,trialvec(:,nsize-ndiv+1:nsize),zm3,.true.)

        ! Build <trialvec|S|trialvec> and <trialvec|trialvec>,
        ! and again we reuse whatever was already available
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

        ! Diagonalise the Ritz matrix and calculate the new estimate for the eigenvectors 
        call diagH2(nsize,blockH,blockS,rd,info)

        if (info.eq.0) then
          call getritzvectors(n_local,nsize,nsize,trialvec,blockH,ritzvec,1,n_local)
          call getritzvectors(n_local,nsize,nsize,Sx,blockH,zm2,1,n_local)
        endif

        deallocate(zm3,BlockS,BlockH)

        ! Calculate the residuals
        res=0d0
        do i=1,ndiv
          resvec(1:npw)=(zm2(1:npw,i)-rd(i)*ritzvec(1:npw,i))
          zsum=zdotc(npw,resvec,1,resvec,1)
          if (dble(zsum).gt.res) res=dble(zsum)
        enddo

        ! Checking the convergence criterion
        keepworking=(abs(res).gt.1.e-8_dp)
        if (mpi_env%is_root .and. (input%groundstate%outputlevel == "high")) &
          write(*,*) sum(rd(1:ndiv)), res, nsize, ii
        evals(1:ndiv)=rd(1:ndiv)
        if ((info /= 0) .or. (.not.keepworking)) ii=nblocks
          
      enddo

      ! Restart the Davidson algorithm
      trialvec=zzero
      trialvec(1:n_local,1:ndiv)=ritzvec(1:n_local,1:ndiv)
      call OverlapX(n_local,npw,ndiv,system,fftmap,zfftcf,trialvec,zm2,.true.)
      Sx(1:n_local,1:ndiv)=zm2(1:n_local,1:ndiv)
      allocate(BlockH(ndiv,ndiv))
      call innerproduct(n_local,ndiv,ndiv,trialvec,zm2,blockH(1:ndiv,1:ndiv))

      extraH=zzero
      extraH(1:ndiv,1:ndiv)=blockH(1:ndiv,1:ndiv)
      extraS=zzero
      do i=1,ndiv
        extraS(i,i)=zone
      enddo
      deallocate(blockH)
    enddo

    ! Store all eigenvectors with the eigenvalue under 0.3. 
    ! The threshold was chosen empirically.
    ! The allocation is conditional because it is done only for the first k-point.

    ! First, determine how many singular components are there FOR THIS k-point
    i = ndiv
    do while ( (i > 1) .and. ( rd(i) > 3.e-1_dp ))
      i = i - 1
    enddo
    if ( i == 0 ) i = 1
    
    ! Second, handle allocations. 
    ! We MUST use the safe upper-bound `ndiv` if allocating for all nkpt.
    if ( .not. ( allocated( singular ) ) ) then
      ! Assuming the SCF scenario with possibly several k-points
      allocate( singular(ngkmax, ndiv, nkpt), source = zzero)
      nsingular = ndiv ! Keep global size at the safe maximum
    elseif (size(singular,3).eq.1) then
      ! 1 kpt? Safe to perfectly fit it to 'i' and reallocate.
      nsingular = i
      deallocate(singular)
      allocate( singular(ngkmax, nsingular, 1), source = zzero )
    endif

    ! Third, copy the result using the LOCAL count 'i'
    ! The rest of the array (up to nsingular) remains perfectly safely as zero.
    singular(1:n_local, 1:i, page) = ritzvec(1:n_local, 1:i)

    if (associated(fftmap%igfft)) deallocate(fftmap%igfft)

    if (mpi_env%is_root .and. (input%groundstate%outputlevel == "high")) then
      write(*,*) calls,'S|trialvec> calls'
      write(*,*) i, 'singular components for this k-point'
      write(*,*) '***********'
      write(*,*) 'highest eigenvalue among singular components', rd(i)
      write(*,*) rd(1:i)
    endif
  endif

  call timesec(tsb)
  if (mpi_env%is_root .and. (input%groundstate%outputlevel == "high")) &
    write(*,*) 'time (singular components):', tsb - tsa
  timefv = timefv + tsb - tsa
end subroutine singularcomponents
