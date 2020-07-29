!BOP
!!ROUTINE: calcpmatgw
!!INTERFACE:
!
subroutine calcpmatgw
!
!!USES:
    use modinput
    use modmain
    use modgw
    use m_getunit
    use modmpi
    use mod_pmat
    use mod_hdf5

!!DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMATVV.OUT} or {\tt PMATCV.OUT}.
!
!!REVISION HISTORY:
!   Created October 2013 (DIN)
!EOP
!BOC
    implicit none
! local variables
    integer    :: ikp, ik, fid, ispn, i, j
    integer(8) :: recl
    real(8)    :: tstart, tend, t0, t1
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: pmv_k(:,:,:), pmc_k(:,:,:)
    complex(8), allocatable :: pmv(:,:,:,:), pmc(:,:,:,:)

    integer    :: k, isym, lspl
    real(8)    :: v(3), v1(3), v2(3), pm(9), sl(3,3), sc(3,3)
    complex(8) :: p(3), o(6)

    integer :: ikstart, ikend
    integer, allocatable :: ikp2rank(:)

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

    !===============================
    ! Initialization
    !===============================
    call init_pmat(input%gw%coreflag=='all', nstdf)

    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfv(nmatmax,nstfv))

    !=========================
    ! Loop over k-points
    !=========================
    allocate(pmv_k(1:nomax,numin:nstdf,3))
    allocate(pmv(1:nomax,numin:nstdf,3,ikstart:ikend))
    pmv(:,:,:,:) = zzero
    if (input%gw%coreflag=='all') then
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
      &                1, nomax, numin, nstdf, pmv_k)
      pmv(:,:,:,ikp) = pmv(:,:,:,ikp) + pmv_k(:,:,:)
      !------------------------------------
      ! core-valence contribution
      !------------------------------------
      if (input%gw%coreflag=='all') then
        call genpmatcv_k(kqset%vkl(:,ik), &
        &                1, ncg, numin, nstdf, pmc_k)
        pmc(:,:,:,ikp) = pmc(:,:,:,ikp) + pmc_k(:,:,:)
      endif

    end do

    deallocate(pmv_k)
    if (input%gw%coreflag=='all') deallocate(pmc_k)
    deallocate(apwalm)
    deallocate(evecfv)
    call clear_pmat()

    !==========================
    ! Write results to files
    !==========================

    ! overwrite existing files
    if (rank==0) then
      call getunit(fid)
      open(fid,File=fname_pmatvv,form='UNFORMATTED',status='REPLACE')
      close(fid)
      if (input%gw%coreflag=='all') then
        call getunit(fid)
        open(fid,File=fname_pmatcv,form='UNFORMATTED',status='REPLACE')
        close(fid)
      end if
    endif
    call barrier

    do ikp = 1, kset%nkpt
      if (rank == ikp2rank(ikp)) then
        call getunit(fid)
        inquire(iolength=recl) pmv(:,:,:,ikp)
        open(fid,File=fname_pmatvv,Action='WRITE',Form='UNFORMATTED',&
             Access='DIRECT',Status='OLD',Recl=recl)
        write(fid,rec=ikp) pmv(:,:,:,ikp)
        close(fid)
        if (input%gw%coreflag=='all') then
          call getunit(fid)
          inquire(iolength=recl) pmc(:,:,:,ikp)
          open(fid,File=fname_pmatcv,Action='WRITE',Form='UNFORMATTED',&
               Access='DIRECT',Status='OLD',Recl=recl)
          write(fid,rec=ikp) pmc(:,:,:,ikp)
          close(fid)
        end if
      end if ! rank
      call barrier
    end do

    deallocate(pmv)
    if (input%gw%coreflag=='all') deallocate(pmc)

    ! timing
    call timesec(tend)
    time_pmat = time_pmat+tend-tstart

    return
end subroutine
!EOC
