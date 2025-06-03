
subroutine calcepsilon(iq,iomstart,iomend)
!
! Compute the RPA dielectric matrix
!
    use modinput
    use modmain, only : zone, zzero, zi
    use modmpi, only: rank
    use modgw
    use modxs,      only : symt2
    use mod_bands, only: nstdf, nomax, numin, eveckpalm, eveckalm, eveck, eveckp
    use precision,  only : i32, dp, long_int
    use iso_c_binding,         only: c_ptr, c_loc, c_f_pointer, c_sizeof, c_size_t
    use mod_device_offload,    only: device_world
    use mod_pointer_remapping, only: remap_fortran_pointer 
    use device_linalg_common_interface, only: zgemm_gpu
    use m_memory_device,       only: allocate_device_memory, deallocate_device_memory, &
                                     bytes_double_complex, bytes_int, get_device_pointer
#include "offload.fpp"
    
    implicit none
    ! input/output
    integer(i32), intent(in) :: iq
    integer(i32), intent(in) :: iomstart, iomend
    ! local
    integer(i32) :: ie1, ie2, ibasis
    integer(i32) :: iom
    integer(i32) :: ik, jk, ispn
    integer(i32) :: im, iop, jop
    integer(i32) :: ndim, mdim, nmdim
    integer(i32) :: nblk, iblk, mstart, mend
    integer(long_int) :: recl
    real(dp)    :: tstart, tend
    real(dp)    :: wto, wlo
    complex(dp) :: head(3,3), f, w
    type(c_ptr) :: minm_cptr
    complex(dp), pointer, contiguous :: minm(:,:,:)
    complex(dp), allocatable :: evecfv(:,:)
    logical :: print_Polarizability
    integer(i32) :: my_device

    call timesec(tstart)

    !=============================
    ! Initialization
    !=============================

    ! Get device id
    my_device = device_world%get_device()

    ! total number of states including the core ones
    if (input%gw%coreflag=='all') then
        ndim = nomax+ncg
    else
        ndim = nomax
    end if
    mdim = nstdf-numin+1

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

    OMP_OFFLOAD target enter data map(alloc: eveckalm, eveckpalm, eveck, eveckp)

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
        
        OMP_OFFLOAD target update to(eveck, eveckp, eveckalm, eveckpalm)


        !=================================================
        ! Loop over m-blocks in M^i_{nm}(\vec{k},\vec{q})
        !=================================================
        do iblk = 1, nblk

            mstart = numin + (iblk-1)*mblksiz
            mend   = min(nstdf, mstart+mblksiz-1)
            nmdim  = ndim * (mend-mstart+1)

            allocate(minmmat(mbsiz,ndim,mstart:mend))
            OMP_OFFLOAD target enter data map(alloc: minmmat)
            msize = sizeof(minmmat)*b2mb

            ! compute M^i_{nm}+M^i_{cm}
            call expand_products(ik, iq, 1, ndim, nomax, mstart, mend, -1, minmmat)

            if (Gamma) then
                ! wings of the dielectric matrix
                OMP_OFFLOAD target update from(minmmat)
                call calcwings(ik, iq, iomstart, iomend, ndim, mstart, mend)
            end if

            ! Body
            call allocate_device_memory(minm_cptr, mbsiz*nmdim*bytes_double_complex, my_device)
            call c_f_pointer(minm_cptr, minm, int([mbsiz,ndim,(mend-mstart+1)],kind=c_size_t))
            ! Remapping the pointer boundaries from Fortran default
            ! TODO(mrm): When supported use lower for c_f_pointer introduced in Fortran 2023
            call remap_fortran_pointer(minm, int([1, 1, mstart], kind=i32), int([mbsiz, ndim, mend], kind=i32))

            do iom = iomstart, iomend
                OMP_OFFLOAD target data map(to: fnm(:,mstart:mend,iom,ik))
                OMP_OFFLOAD target has_device_addr(minm)
                !$omp teams distribute parallel do collapse(3) default(none) private(ie1,ie2,ibasis) &
                !$omp shared(mstart,mend,ndim,mbsiz,minm,fnm,minmmat,iom,ik)
                do ie2 = mstart, mend
                    do ie1 = 1, ndim
                        do ibasis = 1, mbsiz
                            minm(ibasis,ie1,ie2) = fnm(ie1,ie2,iom,ik) * &
                                                   minmmat(ibasis,ie1,ie2)
                        end do
                    end do ! ie1
                end do ! ie2
                !$omp end teams distribute parallel do
                OMP_OFFLOAD end target
                OMP_OFFLOAD end target data
                call zgemm_gpu( 'n', 'c', mbsiz, mbsiz, nmdim, &
                            zone, minm_cptr, mbsiz, get_device_pointer(minmmat,my_device), mbsiz, &
                            zone, get_device_pointer(epsilon(1,1,iom),my_device), mbsiz, device_world)
                call device_world%synchronize()
            end do ! iom
            nullify(minm)
            call deallocate_device_memory(minm_cptr,my_device)
            OMP_OFFLOAD target exit data map(delete: minmmat)
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

    OMP_OFFLOAD target update from(epsilon)

    OMP_OFFLOAD target exit data map(delete: eveck, eveckp, eveckalm, eveckpalm)

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
        if (rank == 0) call writedielt('EPS00', iomend-iomstart+1, freq%freqs(iomstart:iomend), epsh(:,:,iomstart:iomend), 1)
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
            if (rank == 0) call writedielt('EPS00+LAT', iomend-iomstart+1, freq%freqs(iomstart:iomend), epsh(:,:,iomstart:iomend), 1)
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
