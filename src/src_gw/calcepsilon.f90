!BOP
!
!!ROUTINE: calcepsilon
!
!!INTERFACE:
!
subroutine calcepsilon(iq,iomstart,iomend)
!
!!DESCRIPTION:
!
! This subroutine calculates the dielectric matrix in the RPA approximation
! using symmetry
!
!!USES:
    use modinput
    use modmain, only : zone, zzero
    use modgw
    use mod_mpi_gw, only : myrank
    use modxs,      only : symt2
    use m_getunit

!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: iomstart, iomend

!!LOCAL VARIABLES:
    integer(4) :: ia, is, ias
    integer(4) :: ie1, ie2, ie12, ic, icg
    integer(4) :: iom
    integer(4) :: ik, jk, ispn
    integer(4) :: im, jm, iop, jop
    
    integer(4) :: ndim, mdim, nmdim
    
    integer(8) :: recl
      
    real(8)    :: tstart, tend, t0, t1
    
    complex(8) :: coefb
    complex(8) :: head(3,3)
    complex(8), allocatable :: minm(:,:,:)
    complex(8), allocatable :: evecsv(:,:,:)
    
!!EXTERNAL ROUTINES: 
    external zgemm
        
!!REVISION HISTORY:
!
! Created Nov 2013 by DIN
!
!EOP
!BOC
    ! write(*,*)
    ! write(*,*) ' ---- calcepsilon started ----'
    ! write(*,*)
    call timesec(tstart)
    
    ! memory usage info
    msize = sizeof(epsilon)*b2mb
    ! write(*,'(" calcepsilon: rank, size(epsilon) (Mb):",i4,f8.2)') myrank, msize

    !=============================
    ! Initialization
    !=============================
    
    ! constant 'body' prefactor
    coefb = -occmax*zone
    
    ! total number of states including the core ones
    if (input%gw%coreflag=='all') then
      ndim = nomax+ncg
    else
      ndim = nomax
    end if
    mdim = nstdf-numin+1
    nmdim = ndim*mdim
    
    ! arrays to store products of KS eigenvectors with the matching coefficients
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    
    ! memory usage
    msize = (sizeof(eveckalm)+sizeof(eveckpalm)+sizeof(eveck)+sizeof(eveckp))*b2mb
    ! write(*,'(" calcepsilon: rank, size(eigenvectors) (Mb):",i4,f12.2)') myrank, msize
    
    !==================================================
    ! Calculate the q-dependent BZ integration weights
    !==================================================
    call qdepwtet(iq,iomstart,iomend,ndim)

    !==========================
    ! Momentum matrix elements
    !==========================
    if (Gamma) then
      !---------
      ! val-val
      !---------
      if (allocated(pmatvv)) deallocate(pmatvv)
      allocate(pmatvv(nomax,numin:nstdf,3))
      msize = sizeof(pmatvv)*b2mb
      inquire(iolength=recl) pmatvv
      open(fid_pmatvv,File=fname_pmatvv, &
      &    Action='READ',Form='UNFORMATTED',&
      &    Access='DIRECT',Status='OLD',Recl=recl)
      !----------
      ! core-val 
      !----------
      if (input%gw%coreflag=='all') then
        if (allocated(pmatcv)) deallocate(pmatcv)
        allocate(pmatcv(ncg,numin:nstdf,3))
        msize = msize+sizeof(pmatcv)*b2mb
        inquire(iolength=recl) pmatcv
        open(fid_pmatcv,File=fname_pmatcv, &
        &    Action='READ',Form='UNFORMATTED', &
        &    Access='DIRECT',Status='OLD',Recl=recl)
      end if
      ! write(*,'(" calcepsilon: rank, size(pmat) (Mb):",i4,f8.2)') myrank, msize
    end if
    
    allocate(minm(1:mbsiz,1:ndim,numin:nstdf))
    allocate(minmmat(1:mbsiz,1:ndim,numin:nstdf))
    msize = 2*sizeof(minmmat)*b2mb
    ! write(*,'(" calcepsilon: rank, size(minmmat) (Mb):",i4,f12.2)') myrank, msize

    !=================
    ! BZ integration
    !=================
    ! write(*,*) 'calcepsilon: k-point summation'
    do ispn = 1, nspinor
    do ik = 1, kqset%nkpt
    
      ! write(*,*)
      ! write(*,*) 'calcepsilon: rank, (iq, ik):', myrank, iq, ik
    
      ! k-q point
      jk = kqset%kqid(ik,iq)
      
      if (Gamma) then
        !====================================
        ! read the momentum matrix elements
        !====================================
        call getpmatkgw(ik)
        !=================================
        ! Head of the dielectric function
        !=================================
        call calchead(ik,iomstart,iomend,ndim)
      end if
      
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
      call expand_products(ik,iq,1,ndim,nomax,numin,nstdf,-1,minmmat)
        
      if (Gamma) then
        !==============
        ! Wings
        !==============
        call calcwings(ik,iq,iomstart,iomend,ndim,numin,nstdf)
      end if

      !==============
      ! Body
      !==============
      do iom = iomstart, iomend
          
        do ie2 = numin, nstdf
          do ie1 = 1, ndim
            minm(1:mbsiz,ie1,ie2) = fnm(ie1,ie2,iom,ik)* &
            &                       minmmat(1:mbsiz,ie1,ie2)
          end do ! ie1
        end do ! ie2
          
        call zgemm( 'n','c',mbsiz,mbsiz,nmdim, &
        &            coefb,minm,mbsiz,minmmat,mbsiz, &
        &            zone,epsilon(:,:,iom),mbsiz)
          
      end do ! iom
      
    end do ! ik
    end do ! ispn
    
    if (input%gw%selfenergy%secordw) then
      !---------------------------------
      ! Second order screened exchange:
      ! -- store v^{1/2}*P_{0}*v^{1/2}
      !---------------------------------
      vPv(:,:,:) = -epsilon(:,:,:)
      if (Gamma) then
        vPvh(:) = -(epsh(:,1,1)+ epsh(:,2,2)+epsh(:,3,3))/3.d0
        vPvw1(:,:) = -(epsw1(:,:,1)+epsw1(:,:,2)+epsw1(:,:,3))/dsqrt(3.d0)
        vPvw2(:,:) = -(epsw2(:,:,1)+epsw2(:,:,2)+epsw2(:,:,3))/dsqrt(3.d0)
      end if
    end if
    
    !-------------------
    ! Clear memory
    !-------------------
    deallocate(minm)
    deallocate(minmmat)
    if (Gamma) then
      ! deallocate the momentum matrix elements
      deallocate(pmatvv)
      if (input%gw%coreflag=='all') deallocate(pmatcv)
      ! close files
      close(fid_pmatvv)
      if (input%gw%coreflag=='all') close(fid_pmatcv)
    end if
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    
    if (Gamma) then
      !=============================
      ! symmetrize \eps_{00} (head)
      !=============================
      do iom = iomstart, iomend
        head(:,:) = epsh(iom,:,:)
        do iop = 1, 3
          do jop = 1, 3
            call symt2app(iop,jop,1,symt2,head,epsh(iom,iop,jop))
          end do
        end do
      end do ! iom
    end if ! Gamma

    !================================
    ! \epsilon = 1+(-\epsilon) 
    !================================     
    do iom = iomstart, iomend
      if (Gamma) then
        do iop = 1, 3
          epsh(iom,iop,iop) = zone+epsh(iom,iop,iop)
        end do
      end if
      do im = 1, mbsiz
        epsilon(im,im,iom) = zone+epsilon(im,im,iom)
      end do ! im
    end do ! iom
    
    ! timing
    call timesec(tend)
    time_df = time_df+tend-tstart
    
    ! write(*,*) 
    ! write(*,*) ' ---- calcepsilon ended ----'
    ! write(*,*)
      
    return
end subroutine
!EOC
