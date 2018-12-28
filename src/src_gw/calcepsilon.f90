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
    integer(4) :: ibl,nbl,blstart,blend
    
    integer(4) :: ndim, mdim, nmdim
    integer(4) :: nblk, iblk, mstart, mend
    
    integer(8) :: recl
      
    real(8)    :: tstart, tend, t0, t1, ta ,tb
    
    complex(8) :: coefb
    complex(8) :: head(3,3)
    complex(8), allocatable :: minm(:,:,:),mnm(:,:)
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
    
    ! block size
    ! print*, mblksiz
    if (mblksiz >= mdim) then
      nblk = 1
    else
      nblk = mdim / mblksiz
      if ( mod(mdim,mblksiz) /= 0 ) nblk = nblk+1
    end if

    ! arrays to store products of KS eigenvectors with the matching coefficients
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    
    ! memory usage
    ! msize = (sizeof(eveckalm)+sizeof(eveckpalm)+sizeof(eveck)+sizeof(eveckp))*b2mb
    ! write(*,'(" calcepsilon: rank, size(eigenvectors) (Mb):",i4,f12.2)') myrank, msize
    
    !==================================================
    ! Calculate the q-dependent BZ integration weights
    !==================================================
    select case (trim(input%gw%qdepw))
    case('sum')
      call qdepwsum(iq,iomstart,iomend,ndim)
    case('tet')
      call qdepwtet(iq,iomstart,iomend,ndim)
    case default
      stop "Error(calcepsilon): Unknown qdepw method!"
    end select

    !==========================
    ! Momentum matrix elements
    !==========================
    if (Gamma) then
      !---------
      ! val-val
      !---------
      if (allocated(pmatvv)) deallocate(pmatvv)
      allocate(pmatvv(nomax,numin:nstdf,3))
      ! msize = sizeof(pmatvv)*b2mb
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
        ! msize = msize+sizeof(pmatcv)*b2mb
        inquire(iolength=recl) pmatcv
        open(fid_pmatcv,File=fname_pmatcv, &
        &    Action='READ',Form='UNFORMATTED', &
        &    Access='DIRECT',Status='OLD',Recl=recl)
      end if
      ! write(*,'(" calcepsilon: rank, size(pmat) (Mb):",i4,f8.2)') myrank, msize
    end if
    
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
      call getevecsvgw('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
      eveckp = conjg(evecsv(:,:,ispn))
      call getevecsvgw('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
      eveck = evecsv(:,:,ispn)
      deallocate(evecsv)
      call timesec(t1)
      time_io = time_io+t1-t0
        
      call expand_evec(ik,'t')
      call expand_evec(jk,'c')

      !=================================================
      ! Loop over m-blocks in M^i_{nm}(\vec{k},\vec{q})
      !=================================================
      do iblk = 1, nblk

        ! call timesec(ta)      
        mstart = numin+(iblk-1)*mblksiz
        mend   = min(nstdf,mstart+mblksiz-1)
        nmdim  = ndim*(mend-mstart+1)
        ! print*, iblk, nblk, mstart, mend

        allocate(minmmat(mbsiz,ndim,mstart:mend))
        ! msize = sizeof(minmmat)*b2mb
        ! write(*,'(" calcepsilon: rank, size(minmmat) (Mb):",3i4,f12.2)') myrank, mstart, mend, msize

        ! Calculate M^i_{nm}+M^i_{cm}
        call expand_products(ik, iq, 1, ndim, nomax, mstart, mend, -1, minmmat)
        
        if (Gamma) then
          !==============
          ! Wings
          !==============
          call calcwings(ik, iq, iomstart, iomend, ndim, mstart, mend)
        end if

        !==============
        ! Body
        !==============
        allocate(minm(mbsiz,ndim,mstart:mend))
        do iom = iomstart, iomend

          do ie2 = mstart, mend
            do ie1 = 1, ndim
              minm(1:mbsiz,ie1,ie2) = fnm(ie1,ie2,iom,ik) * &
              &                       minmmat(1:mbsiz,ie1,ie2)
            end do ! ie1
          end do ! ie2

          call zgemm( 'n', 'c', mbsiz, mbsiz, nmdim, &
          &            coefb, minm, mbsiz, minmmat, mbsiz, &
          &            zone, epsilon(:,:,iom), mbsiz)

        end do ! iom

        deallocate(minm, minmmat)
        !call timesec(tb)      
        !write(*,*) iblk,' time:',tb-ta

      end do ! iblk
      
    end do ! ik
    end do ! ispn
    
    !-------------------
    ! Clear memory
    !-------------------
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
            call symt2app(iop, jop, 1, symt2, head, epsh(iom,iop,jop))
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
