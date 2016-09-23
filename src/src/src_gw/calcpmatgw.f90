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
    integer    :: ik, ikp, fid, ispn, i
    integer(8) :: recl
    real(8) :: tstart, tend, t0, t1
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecsv(:,:,:)
    complex(8), allocatable :: pmv_k(:,:,:), pmc_k(:,:,:)
    complex(8), allocatable :: pmv(:,:,:,:), pmc(:,:,:,:)
    
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
    do ik = ikstart, ikend
      ikp2rank(ik) = rank
    end do ! ik

    !===============================
    ! Initialization
    !===============================
    call init_pmat(input%gw%coreflag=='all')

    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecsv(nmatmax,nstsv,nspinor))
    
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
      call match(Gkset%ngk(1,ik),Gkset%gkc(:,1,ik), &
      &          Gkset%tpgkc(:,:,1,ik),Gkset%sfacgk(:,:,1,ik), &
      &          apwalm)
      
      !-------------------------------------------
      ! get the eigenvectors and values from file
      !-------------------------------------------
      call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
      
      do ispn = 1, nspinor
        call genevecalm(Gkset%ngk(1,ik),nstsv,evecsv(:,:,ispn),apwalm)
        !------------------------------------      
        ! valence-valence contribution
        !------------------------------------
        call genpmatvv_k(Gkset%ngk(1,ik), Gkset%igkig(:,1,ik), &
        &                Gkset%vgkc(:,:,1,ik), nstsv, evecsv(:,:,ispn), &
        &                1, nomax, numin, nstdf, pmv_k)
        pmv(:,:,:,ikp) = pmv(:,:,:,ikp)+pmv_k(:,:,:)
        !------------------------------------      
        ! core-valence contribution
        !------------------------------------
        if (input%gw%coreflag=='all') then
          call genpmatcv_k(kqset%vkl(:,ik), &
          &                1, ncg, numin, nstdf, pmc_k)
          pmc(:,:,:,ikp) = pmc(:,:,:,ikp)+pmc_k(:,:,:)
        endif
      end do ! ispn

    end do ! ikp
    
    deallocate(pmv_k)
    if (input%gw%coreflag=='all') deallocate(pmc_k)
    deallocate(apwalm)
    deallocate(evecsv)
    call clear_pmat
    
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
      if (rank==ikp2rank(ikp)) then
#ifdef _HDF5_
        write(cik,'(I4.4)') ikp
        path = "/kpoints/"//trim(adjustl(cik))
        if (.not.hdf5_exist_group(fgwh5,"/kpoints",cik)) &
        &  call hdf5_create_group(fgwh5,"/kpoints",cik)
        call hdf5_write(fgwh5,path,"pmatvv", &
        &               pmv(1,numin,1,ikp),(/nomax,nstdf-numin+1,3/))
        if (input%gw%coreflag=='all') then
          call hdf5_write(fgwh5,path,"pmatcv", &
          &               pmc(1,numin,1,ikp),(/ncg,nstdf-numin+1,3/))
        end if
#endif
        call getunit(fid)
        inquire(iolength=recl) pmv(:,:,:,ikp)
        open(fid,File=fname_pmatvv,Action='WRITE',Form='UNFORMATTED',&
        &    Access='DIRECT',Status='OLD',Recl=recl)
        write(fid,rec=ikp) pmv(:,:,:,ikp)
        close(fid)
        if (input%gw%coreflag=='all') then
          call getunit(fid)
          inquire(iolength=recl) pmc(:,:,:,ikp)
          open(fid,File=fname_pmatcv,Action='WRITE',Form='UNFORMATTED',&
          &    Access='DIRECT',Status='OLD',Recl=recl)
          write(fid,rec=ikp) pmc(:,:,:,ikp)
          close(fid)
        end if
      end if ! rank
      call barrier
    end do ! ikp

    deallocate(pmv)
    if (input%gw%coreflag=='all') deallocate(pmc)
    
    ! timing
    call timesec(tend)
    time_pmat = time_pmat+tend-tstart
    
    return
end subroutine
!EOC
