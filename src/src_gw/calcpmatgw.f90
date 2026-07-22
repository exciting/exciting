!BOP
!!ROUTINE: calcpmatgw
!!INTERFACE:
!
subroutine calcpmatgw
!
!!USES:
    use modinput
    use constants, only: zzero
    use mod_APW_LO, only: apwordmax
    use mod_muffin_tin, only: lmmaxapw
    use mod_Gkvector, only: ngkmax
    use mod_atoms, only: natmtot
    use mod_eigensystem, only: nmatmax
    use mod_core_states,   only : ncg
    use mod_large_io, only: inquire_large, open_direct_unformatted_large
    use mod_pmat, only: init_pmat, genevecalm, genpmatvv_k, genpmatcv_k, clear_pmat
    use mod_dielectric_function, only: fname_pmatvv, fname_pmatcv
    use modgw, only: kset, kqset, Gkqset, time_pmat 
    use mod_bands, only: nomax, numin, nstdf
    use mod_eigenvalue_occupancy, only: nstfv
    use modmpi, only: rank, firstofset, lastofset, barrier 

!!DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMATVV.OUT} or {\tt PMATCV.OUT}.
!
!!REVISION HISTORY:
!   Created October 2013 (DIN)
!EOP
!BOC
    use precision, only: i32, long_int, dp
    implicit none
! local variables
    integer(i32)  :: ikp, ik, fid, ispn, i, j
    integer(long_int) :: recl
    real(dp)    :: tstart, tend, t0, t1
    complex(dp), allocatable :: apwalm(:,:,:,:)
    complex(dp), allocatable :: evecfv(:,:)
    complex(dp), allocatable :: pmv_k(:,:,:), pmc_k(:,:,:)
    complex(dp), allocatable :: pmv(:,:,:,:), pmc(:,:,:,:)

    integer(i32) :: k, isym, lspl
    real(dp)     :: v(3), v1(3), v2(3), pm(9), sl(3,3), sc(3,3)
    complex(dp)  :: p(3), o(6)

    integer(i32) :: ikstart, ikend
    integer(i32), allocatable :: ikp2rank(:)

    logical(i32) :: coreflag_is_all

    call timesec(tstart)

#ifdef MPI
    ikstart = firstofset(rank,kset%nkpt)
    ikend = lastofset(rank,kset%nkpt)
#else
    ikstart = 1
    ikend = kset%nkpt
#endif

    allocate(ikp2rank(kset%nkpt))
    ikp2rank = -1
    do ikp = ikstart, ikend
      ikp2rank(ikp) = rank
    end do ! ik

    coreflag_is_all = input%gw%coreflag=='all'

    !===============================
    ! Initialization
    !===============================
    call init_pmat(coreflag_is_all, nstdf)

    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfv(nmatmax,nstfv))

    !=========================
    ! Loop over k-points
    !=========================
    allocate(pmv_k(1:nomax,numin:nstdf,3))
    allocate(pmv(1:nomax,numin:nstdf,3,ikstart:ikend))
    pmv(:,:,:,:) = zzero
    if (coreflag_is_all) then
      allocate(pmc_k(1:ncg,numin:nstdf,3))
      allocate(pmc(1:ncg,numin:nstdf,3,ikstart:ikend))
      pmc(:,:,:,:) = zzero
    end if

    do ikp = ikstart, ikend
      ik = kset%ikp2ik(ikp)
      !-------------------------------------------
      ! find the matching coefficients
      !-------------------------------------------
      call match(Gkqset%ngk(1,ik), Gkqset%gkc(:,1,ik), &
      &          Gkqset%tpgkc(:,:,1,ik), Gkqset%sfacgk(:,:,1,ik), &
      &          apwalm)
      !-------------------------------------------
      ! get the eigenvectors and values from file
      !-------------------------------------------
      call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
      call genevecalm(Gkqset%ngk(1,ik), nstdf, evecfv(:,1:nstdf), apwalm)
      !------------------------------------
      ! valence-valence contribution
      !------------------------------------
      call genpmatvv_k(Gkqset%ngk(1,ik), Gkqset%igkig(:,1,ik), &
      &                Gkqset%vgkc(:,:,1,ik), nstdf, evecfv(:,1:nstdf), &
      &                1_i32, nomax, numin, nstdf, pmv_k)
      pmv(:,:,:,ikp) = pmv(:,:,:,ikp) + pmv_k(:,:,:)
      !------------------------------------
      ! core-valence contribution
      !------------------------------------
      if (coreflag_is_all) then
        call genpmatcv_k(kqset%vkl(:,ik), &
        &                1_i32, ncg, numin, nstdf, pmc_k)
        pmc(:,:,:,ikp) = pmc(:,:,:,ikp) + pmc_k(:,:,:)
      endif

    end do

    deallocate(pmv_k)
    if (coreflag_is_all) deallocate(pmc_k)
    deallocate(apwalm)
    deallocate(evecfv)
    call clear_pmat()

    !==========================
    ! Write results to files
    !==========================

    ! overwrite existing files
    if (rank==0) then
      open(newunit=fid,File=fname_pmatvv,form='UNFORMATTED',status='REPLACE')
      close(fid)
      if (coreflag_is_all) then
        open(newunit=fid,File=fname_pmatcv,form='UNFORMATTED',status='REPLACE')
        close(fid)
      end if
    endif
    call barrier

    do ikp = 1, kset%nkpt
      if (rank == ikp2rank(ikp)) then
        call inquire_large( recl, pmv(:,:,:,ikp) )
        call open_direct_unformatted_large( fid, fname_pmatvv, "write", recl, "old" )
        write(fid,rec=ikp) pmv(:,:,:,ikp)
        close(fid)
        if (coreflag_is_all) then
          call inquire_large( recl, pmc(:,:,:,ikp) )
          call open_direct_unformatted_large( fid, fname_pmatcv, "write", recl, "old" )
          write(fid,rec=ikp) pmc(:,:,:,ikp)
          close(fid)
        end if
      end if ! rank
      call barrier
    end do

    deallocate(pmv)
    if (coreflag_is_all) deallocate(pmc)

    ! timing
    call timesec(tend)
    time_pmat = time_pmat+tend-tstart

    return
end subroutine
!EOC

