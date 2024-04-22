
subroutine calcepsilon(iq,iomstart,iomend)
!
! Compute the RPA dielectric matrix
!
    use modinput
    use modmain, only : zone, zzero, zi
    use modgw
    use mod_mpi_gw, only : myrank
    use modxs,      only : symt2
    
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
    logical :: print_Polarizability

    external zgemm

    call timesec(tstart)

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
    
    print_Polarizability = associated(input%gw%taskGroup)
    if( print_Polarizability ) print_Polarizability = associated(input%gw%taskGroup%epsilon)
    if( print_Polarizability ) print_Polarizability = input%gw%taskGroup%epsilon%printPolarizabilityFactor
    if( print_Polarizability ) call write_fnm_to_file( iq, fnm, lbound(fnm), input%gw%taskGroup%outputFormat )

    !=================
    ! BZ integration
    !=================
    do ik = 1, kqset%nkpt

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

            mstart = numin + (iblk-1)*mblksiz
            mend   = min(nstdf, mstart+mblksiz-1)
            nmdim  = ndim * (mend-mstart+1)

            allocate(minmmat(mbsiz,ndim,mstart:mend))
            msize = sizeof(minmmat)*b2mb

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

        end do ! iblk

    end do ! ik

    !-------------------
    ! Clear memory
    !-------------------
    if (Gamma) then
        ! deallocate the momentum matrix elements
        deallocate(pmatvv)
        if (input%gw%coreflag=='all') deallocate(pmatcv)
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

contains 
!> (private) Write the polarizability factor into an output file
subroutine write_fnm_to_file( idx_qpoint, polarizability_factor, lbounds, file_format )
    use gw_io, only: build_file_name, write_to_file
    use precision, only: dp, i32, str_64

    implicit none

    !> q-point index
    integer(i32), intent(in)  :: idx_qpoint 
    !> lower bounds of `polarizability_factor`
    integer(i32), intent(in)  :: lbounds(4)
    !> Polarizability factor
    complex(dp), intent(in) :: polarizability_factor(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):)
    !> Format of output file
    character(len=*), intent(in) :: file_format

    character(len=*), parameter :: file_name_polarizability_factor = 'POLARIZABILITY_FACTOR_Q'
  
    character(len=str_64) :: file_name
  
    call build_file_name( file_name_polarizability_factor, idx_qpoint, file_name )
    call write_to_file( file_name, polarizability_factor, lbounds, file_format )
  
end subroutine  
end subroutine
