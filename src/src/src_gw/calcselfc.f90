!BOP
!!ROUTINE: calcselfc
!
!!INTERFACE: 
!
subroutine calcselfc(iq)
!      
!!DESCRIPTION:
!
! This subroutine calculates the q-dependent correaltion self-energy
!
!!USES:
    use modinput
    use modmain,    only : nstsv, apwordmax, lmmaxapw, natmtot, nspnfv, &
    &                      zzero, nmatmax
    use modgw
    use mod_mpi_gw, only : myrank, myrank_col
    use m_getunit
    
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq
    
!!LOCAL VARIABLES:            
    integer(4) :: ik, ikp, jk, ispn
    integer(4) :: mdim
    integer(4) :: fid
    character(120) :: fname_mwm
    real(8) :: tstart, tend, t0, t1
    complex(8), allocatable :: evecsv(:,:,:)
    character(len=10), external :: int2str

!!REVISION HISTORY:
!
! Created Nov 2013 by DIN
!
!EOP
!BOC
    ! if (myrank==0) then
    !   write(*,*)
    !   write(*,*) ' ---- calcselfc started ----'
    !   write(*,*)
    ! end if
    call timesec(tstart)
    
    !------------------------
    ! total number of states
    !------------------------
    if (input%gw%coreflag=='all') then
      mdim = nstse+ncg
    else
      mdim = nstse
    end if
    
    !-------------------------------------------
    ! products M*W^c*M
    !-------------------------------------------
    if (myrank_col==0) then
      allocate(mwm(ibgw:nbgw,1:mdim,1:freq%nomeg))
      msize = sizeof(mwm)*b2mb
      ! write(*,'(" calcselfc: size(mwm) (Mb):",f12.2)') msize
      !----------------------------
      ! q-dependent M*W*M products
      !----------------------------
      if (input%gw%taskname.eq.'gw0') then
        fname_mwm = 'MWM'//'-q'//trim(int2str(iq))//'.OUT'
        call getunit(fid)
        open(fid,File=fname_mwm,Action='Write',Form='Unformatted')
      end if ! myrank_col
    end if
    
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    allocate(minmmat(1:mbsiz,ibgw:nbgw,1:mdim))
    msize = sizeof(minmmat)*b2mb
    ! write(*,'(" calcepsilon: rank, size(minmmat) (Mb):",i4,f12.2)') myrank, msize
    
    !================================
    ! loop over irreducible k-points
    !================================
    do ispn = 1, nspinor
    do ikp = 1, kset%nkpt
    
      ! write(*,*)
      ! write(*,*) 'calcselfc: k-point loop ikp=', ikp
    
      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector 
      jk = kqset%kqid(ik,iq)
      
      ! get KS eigenvectors
      call timesec(t0)
      allocate(evecsv(nmatmax,nstsv,nspinor))
      call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
      eveckp = conjg(evecsv(:,:,ispn))
      call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
      eveck = evecsv(:,:,ispn)
      deallocate(evecsv)
      call timesec(t1)
      time_io = time_io+t1-t0
        
      ! Calculate M^i_{nm}+M^i_{cm}
      call expand_evec(ik,'t')
      call expand_evec(jk,'c')
      call expand_products(ik,iq,ibgw,nbgw,-1,1,mdim,nstse,minmmat)

      !================================================================
      ! Calculate weight(q)*Sum_ij{M^i*W^c_{ij}(k,q;\omega)*conjg(M^j)}
      !================================================================
      call calcmwm(ibgw,nbgw,1,mdim)
      
      !===========================================================
      ! Calculate the contribution to the correlation self-energy
      !===========================================================
      if (myrank_col==0) then
        
        if (input%gw%taskname=='cohsex') then
          call calcselfc_cohsex(ikp,iq,mdim)
        else
          call calcselfc_freqconv(ikp,iq,mdim)
        end if
          
        if (input%gw%taskname=='gw0') then
          ! store M*W*M in files
          call timesec(t0)
          write(fid) mwm
          call timesec(t1)
          time_io = time_io+t1-t0
        end if

      end if
      
      !if (input%gw%selfenergy%secordw) then 
      !---------------------------------
      ! Second order screened exchange
      !---------------------------------
      !  call calcsecordwselfc(ikp,iq,mdim)
      !end if
    
    end do ! ikp
    end do ! ispn
    
    deallocate(minmmat)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    
    ! delete MWM
    if (myrank_col==0) then
      deallocate(mwm)
      ! and close the file
      if (input%gw%taskname.eq.'gw0') close(fid)
    end if
    
    ! timing
    call timesec(tend)
    time_selfc = time_selfc+tend-tstart
    
    ! if (myrank==0) then
    !   write(*,*)
    !   write(*,*) ' ---- calcselfc ended ----'
    !   write(*,*)
    ! end if
    
    return
end subroutine
                        
