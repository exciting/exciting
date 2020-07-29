
subroutine calcepsilon(iq,iomstart,iomend)
!
! Compute the RPA dielectric matrix
!
    use modinput
    use modmain, only : zone, zzero, zi
    use modgw
    use mod_mpi_gw, only : myrank
    use modxs,      only : symt2
    use m_getunit
    implicit none
    ! input/output
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: iomstart, iomend
    ! local
    integer(4) :: ie1, ie2
    integer(4) :: iom
    integer(4) :: ik, jk, ispn
    integer(4) :: im, iop, jop
    integer(4) :: ndim, mdim, nmdim
    integer(4) :: nblk, iblk, mstart, mend
    integer(8) :: recl
    real(8)    :: tstart, tend
    real(8)    :: wto, wlo
    complex(8) :: head(3,3), f, w
    complex(8), allocatable :: minm(:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    external zgemm

    call timesec(tstart)

    ! memory usage info
    ! msize = sizeof(epsilon)*b2mb
    ! write(*,'(" calcepsilon: rank, size(epsilon) (Mb):",i4,f8.2)') myrank, msize

    !=============================
    ! Initialization
    !=============================

    ! total number of states including the core ones
    if (input%gw%coreflag=='all') then
        ndim = nomax+ncg
    else
        ndim = nomax
    end if
    mdim = nstdf-numin+1
    nmdim = ndim*mdim

    ! block size
    if (mblksiz >= mdim) then
      nblk = 1
    else
      nblk = mdim / mblksiz
      if ( mod(mdim,mblksiz) /= 0 ) nblk = nblk+1
    end if

    ! arrays to store products of KS eigenvectors with the matching coefficients
    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstfv))
    allocate(eveckp(nmatmax,nstfv))
    ! msize = (sizeof(eveckalm)+sizeof(eveckpalm)+sizeof(eveck)+sizeof(eveckp))*b2mb
    ! write(*,'(" calcepsilon: rank, size(eigenvectors) (Mb):",i4,f12.2)') myrank, msize

    !==================================================
    ! Calculate the q-dependent BZ integration weights
    !==================================================
    select case (trim(input%gw%qdepw))
    case('sum')
        call qdepwsum(iq, iomstart, iomend, ndim)
    case('tet')
        call qdepwtet(iq, iomstart, iomend, ndim)
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
            inquire(iolength=recl) pmatcv
            open(fid_pmatcv,File=fname_pmatcv, &
            &    Action='READ',Form='UNFORMATTED', &
            &    Access='DIRECT',Status='OLD',Recl=recl)
        end if
    end if

    !=================
    ! BZ integration
    !=================
    ! write(*,*)
    do ik = 1, kqset%nkpt
        ! write(*,*) 'calcepsilon: rank, (iq, ik):', myrank, iq, ik

        ! k-q point
        jk = kqset%kqid(ik, iq)

        if (Gamma) then
            ! read the momentum matrix elements
            call getpmatkgw(ik)
            ! and compute the head of the dielectric function
            call calchead(ik, iomstart, iomend, ndim)
        end if

        ! get KS eigenvectors
        allocate(evecfv(nmatmax,nstfv))
        call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
        eveckp = conjg(evecfv)
        call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
        eveck = evecfv
        deallocate(evecfv)

        ! compute products \sum_G C_{k}n * A_{lm}
        call expand_evec(ik,'t')
        call expand_evec(jk,'c')

        !=================================================
        ! Loop over m-blocks in M^i_{nm}(\vec{k},\vec{q})
        !=================================================
        do iblk = 1, nblk

            ! call timesec(ta)
            mstart = numin + (iblk-1)*mblksiz
            mend   = min(nstdf, mstart+mblksiz-1)
            nmdim  = ndim * (mend-mstart+1)
            ! print*, iblk, nblk, mstart, mend

            allocate(minmmat(mbsiz,ndim,mstart:mend))
            msize = sizeof(minmmat)*b2mb
            ! write(*,'(" calcepsilon: rank, size(minmmat) (Mb):",3i4,f12.2)') myrank, mstart, mend, msize

            ! compute M^i_{nm}+M^i_{cm}
            call expand_products(ik, iq, 1, ndim, nomax, mstart, mend, -1, minmmat)

            if (Gamma) then
                ! wings of the dielectric matrix
                call calcwings(ik, iq, iomstart, iomend, ndim, mstart, mend)
            end if

            ! Body
            allocate(minm(mbsiz,ndim,mstart:mend))
            do iom = iomstart, iomend
                do ie2 = mstart, mend
                    do ie1 = 1, ndim
                        minm(1:mbsiz,ie1,ie2) = fnm(ie1,ie2,iom,ik) * &
                                                minmmat(1:mbsiz,ie1,ie2)
                    end do ! ie1
                end do ! ie2
                call zgemm( 'n', 'c', mbsiz, mbsiz, nmdim, &
                            zone, minm, mbsiz, minmmat, mbsiz, &
                            zone, epsilon(:,:,iom), mbsiz)
            end do ! iom
            deallocate(minm)

            deallocate(minmmat)
            ! call timesec(tb)
            ! write(*,*) iblk,' time:',tb-ta

        end do ! iblk

    end do ! ik

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
        ! symmetrize \eps_{00} (head)
        do iom = iomstart, iomend
            head(:,:) = epsh(:,:,iom)
            do iop = 1, 3
            do jop = 1, 3
                call symt2app(iop, jop, 1, symt2, head, epsh(iop,jop,iom))
            end do
            end do
        end do ! iom
    end if ! Gamma

    !=================================
    ! \epsilon = \delta_ij - \epsilon
    !=================================
    do iom = iomstart, iomend
        if (Gamma) then
            do iop = 1, 3
                epsh(iop,iop,iom) = zone - epsh(iop,iop,iom)
            end do
        end if
        epsilon(:,:,iom) = -epsilon(:,:,iom)
        do im = 1, mbsiz
            epsilon(im,im,iom) = zone + epsilon(im,im,iom)
        end do ! im
    end do ! iom

    ! Compute contributions due to polar phonons
    if (Gamma) then
        if (myrank == 0) call writedielt('EPS00', iomend-iomstart+1, freq%freqs(iomstart:iomend), epsh(:,:,iomstart:iomend), 1)
        if (input%gw%eph == 'polar') then
            ! call eph_polar(iomend-iomstart+1, cmplx(0.d0,freq%freqs(iomstart:iomend),8), epsh(:,:,iomstart:iomend))
            wlo = input%gw%wlo
            wto = input%gw%wto
            do iom = iomstart, iomend
                w = zi*freq%freqs(iom)
                ! f = (wlo**2-w**2) / (wto**2-w**2)
                f = wlo**2 / wto**2
                epsh(:,:,iom) = epsh(:,:,iom)*f
            end do
            if (myrank == 0) call writedielt('EPS00+LAT', iomend-iomstart+1, freq%freqs(iomstart:iomend), epsh(:,:,iomstart:iomend), 1)
        end if
    end if

    ! timing
    call timesec(tend)
    time_df = time_df+tend-tstart

    return
end subroutine
