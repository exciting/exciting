!> This module contains all the procedures related with the wings and heads of the dielectric/polarizability functions
module mod_head_and_wings

    use constants, only: fourpi, zzero, zone, pi
    use precision, only: i32, dp
    use mod_misc_gw,           only: vi
#include "offload.fpp"
    private
    public :: calchead, calcwings

contains

    !> Computes the head.
    !> Note: this procedure does not reset the head to zero.
    !> It accumulates contributions over all k-points.
    subroutine calchead(ik, first_empty_band, last_empty_band, iomstart, iomend, ndim, head)

        use modinput,  only: input
        use modmain,   only:  evalcr, idxas
        use mod_bands, only: numin, nstdf, nomax, evalfv, metallic
        use modgw, only: fdebug, kset, fnm, kwfer, time_dfhead, freq
        use mod_core_states,   only: corind, ncg
        use mod_dielectric_function, only: pmatvv, pmatcv 
        implicit none

        ! input/output
        integer(i32), intent(in) :: ik
        !> Index of the first unoccupied state
        integer(i32), intent(in) :: first_empty_band
        !> Index of the last unoccupied state
        integer(i32), intent(in) :: last_empty_band
        integer(i32), intent(in) :: iomstart, iomend
        integer(i32), intent(in) :: ndim
        complex(dp), intent(inout) :: head(3,3,iomstart:iomend)

        ! local
        integer(i32) :: ic, icg
        integer(i32) :: ia, is, ias
        integer(i32) :: ie1, ie2
        integer(i32) :: ikp, iop, jop
        integer(i32) :: iom
        real(dp) :: edif    ! energy difference
        real(dp) :: tstart, tend
        complex(dp) :: coefh
        complex(dp) :: pnm, zsum
        real(dp), parameter :: tol = 1.e-6_dp

        if (input%gw%debug) write(fdebug,*) ' ---- calchead started ----'
        call timesec(tstart)

        ! position in the non-reducied grid
        ikp = kset%ik2ikp(ik)

        ! constant prefactor
        coefh = cmplx(4.0_dp*pi*vi, 0.0_dp, kind=dp)

        ! loop over tensor components
        do jop = 1, 3
        do iop = 1, 3

            ! loop over frequencies
            do iom = iomstart, iomend

                ! Inter-band contribution
                zsum = zzero
                do ie2 = first_empty_band, last_empty_band
                  do ie1 = 1, ndim
                    if (ie1 <= nomax) then
                        ! valence-valence
                        edif = evalfv(ie2,ikp)-evalfv(ie1,ikp)
                        if (abs(edif) > 1.0e-6_dp) then
                            pnm = pmatvv(ie1,ie2,iop)*conjg(pmatvv(ie1,ie2,jop))
                            zsum = zsum + fnm(ie1,ie2,iom)*pnm/(edif*edif)
                        end if
                    else
                        ! core-valence
                        icg  = ie1 - nomax
                        is   = corind(icg,1)
                        ia   = corind(icg,2)
                        ias  = idxas(ia,is)
                        ic   = corind(icg,6)
                        edif = evalfv(ie2,ikp) - evalcr(ic,ias)
                        if (abs(edif) > tol) then
                            pnm = pmatcv(icg,ie2,iop)*conjg(pmatcv(icg,ie2,jop))
                            zsum = zsum + fnm(ie1,ie2,iom)*pnm/(edif*edif)
                        end if
                    end if
                  end do ! ie2
                end do ! ie1
                head(iop,jop,iom) = head(iop,jop,iom) + coefh*zsum

                !-------------------------
                ! Intra-band contribution
                !-------------------------
                if (metallic) then
                    zsum = zzero
                    do ie1 = first_empty_band, min(nomax, last_empty_band)
                        pnm = pmatvv(ie1,ie1,iop)*conjg(pmatvv(ie1,ie1,jop))
                        zsum = zsum + kwfer(ie1,ik)*pnm
                    enddo
                    ! for imaginary frequency, a negative sign is needed
                    if (trim(freq%fconv) == 'imfreq') zsum = -zsum
                    head(iop,jop,iom) = head(iop,jop,iom) + &
                                        coefh*zsum/(freq%freqs(iom)**2)
                end if ! metallic

            end do ! iom

        end do ! iop
        end do ! jop

        ! timing
        call timesec(tend)
        time_dfhead = time_dfhead+tend-tstart
        if (input%gw%debug) then
            write(fdebug,*) ' ---- calchead finished ----'
            write(fdebug,*) ' ik = ', ik
            do iom = iomstart, iomend
                write(fdebug,*) iom, head(1,1,iom)
            end do
        end if

        return
    end subroutine calchead


    !> Computes the wings
    !> Note: this procedure does not reset the head to zero.
    !> It accumulates contributions over all k-points.
    subroutine calcwings(ik, iq, iomstart, iomend, ndim, mstart, mend, mbsiz, minmmat, wing1, wing2)

        use modmain,   only : evalcr, idxas
        use mod_bands, only: nomax, evalfv
        use modgw, only: kset, fnm, freq, time_dfwing
        use mod_core_states,   only: corind, ncg
        use mod_dielectric_function, only: pmatvv, pmatcv
        use iso_c_binding,         only: c_ptr, c_loc, c_f_pointer, c_sizeof, c_size_t
        use mod_device_offload,    only: device_world
        use mod_pointer_remapping, only: remap_fortran_pointer
        use device_linalg_common_interface, only: zgemm_gpu
        use m_memory_device,       only: allocate_device_memory, deallocate_device_memory, &
                                     bytes_double_complex, get_device_pointer

        implicit none

        ! input variables
        integer(i32), intent(in) :: ik
        integer(i32), intent(in) :: iq
        integer(i32), intent(in) :: iomstart, iomend
        integer(i32), intent(in) :: ndim
        integer(i32), intent(in) :: mstart, mend
        integer(i32), intent(in) :: mbsiz
        complex(dp), intent(in)  :: minmmat(mbsiz,ndim,mstart:mend)
        complex(dp), intent(inout) :: wing1(mbsiz,3,iomstart:iomend)
        complex(dp), intent(inout) :: wing2(mbsiz,3,iomstart:iomend)

        integer(i32) :: nmdim
        integer(i32) :: ie1, ie2, imix
        integer(i32) :: is, ia, ic, icg, ias
        integer(i32) :: iop, iom, ikp
        integer(i32) :: my_device
        real(dp) :: edif
        complex(dp) :: coefw
        real(dp) :: t0, t1

        type(c_ptr) :: pm1_cptr, pm2_cptr, tmat1_cptr, tmat2_cptr
        complex(dp), pointer, contiguous :: pm1(:,:,:), pm2(:,:,:), tmat1(:,:,:), tmat2(:,:,:)

        ! Timer
        call timesec(t0)

        ! Get my device id
        my_device = device_world%get_device()

        ! wings prefactors
        coefw = cmplx(-sqrt(fourpi*vi), 0.0_dp, kind=dp)

        ! irreducible k-point index
        ikp = kset%ik2ikp(ik)

        nmdim = ndim*(mend-mstart+1)

        call allocate_device_memory(pm1_cptr, 3*nmdim*bytes_double_complex, my_device)
        call allocate_device_memory(pm2_cptr, 3*nmdim*bytes_double_complex, my_device)
        call c_f_pointer(pm1_cptr, pm1, int([ndim,(mend-mstart+1), 3],kind=c_size_t))
        call c_f_pointer(pm2_cptr, pm2, int([ndim,(mend-mstart+1), 3],kind=c_size_t))
        ! Remapping the pointer boundaries from Fortran default
        ! TODO(mrm): When supported use lower for c_f_pointer introduced in Fortran 2023
        call remap_fortran_pointer(pm1, int([1, mstart, 1], kind=i32), int([ndim, mend, 3], kind=i32))
        call remap_fortran_pointer(pm2, int([1, mstart, 1], kind=i32), int([ndim, mend, 3], kind=i32))

        ! Iterate over the Cartesian coordinates
        ! and prepare the momentum matrix elements
        OMP_OFFLOAD target has_device_addr(pm1, pm2) map(to: pmatvv, pmatcv)
        !$omp teams distribute parallel do collapse(3) &
        !$omp default(none) private(iop, ie2, ie1, edif, icg, is, ia, ic, ias) &
        !$omp shared(evalfv, ikp, pmatvv, pmatcv, nomax, evalcr, pm1, pm2, corind, mstart, mend, ndim, idxas)
        do iop = 1, 3
            do ie2 = mstart, mend
                do ie1 = 1, ndim
                    ! Valence states
                    if (ie1 <= nomax) then
                        edif = evalfv(ie1,ikp) - evalfv(ie2,ikp)
                        pm1(ie1,ie2,iop) = merge(pmatvv(ie1,ie2,iop) / edif, zzero, abs(edif) > 1.0e-6_dp)
                        pm2(ie1,ie2,iop) = conjg(pm1(ie1,ie2,iop))
                    ! Core states
                    else
                        icg = ie1 - nomax
                        is  = corind(icg,1)
                        ia  = corind(icg,2)
                        ic  = corind(icg,6)
                        ias = idxas(ia,is)
                        edif = evalcr(ic,ias) - evalfv(ie2,ikp)
                        pm1(ie1,ie2,iop) = merge(pmatcv(icg,ie2,iop)/edif, zzero, abs(edif) > 1.0e-6_dp)
                        pm2(ie1,ie2,iop) = conjg(pm1(ie1,ie2,iop))
                    end if
                end do
            end do
        end do
        !$omp end teams distribute parallel do
        OMP_OFFLOAD end target

        ! Build the expansion coefficient multiplied by the occupation and frequency factors
        call allocate_device_memory(tmat1_cptr, mbsiz*nmdim*bytes_double_complex, my_device)
        call allocate_device_memory(tmat2_cptr, mbsiz*nmdim*bytes_double_complex, my_device)
        call c_f_pointer(tmat1_cptr, tmat1, int([mbsiz,ndim,(mend-mstart+1)],kind=c_size_t))
        call c_f_pointer(tmat2_cptr, tmat2, int([mbsiz,ndim,(mend-mstart+1)],kind=c_size_t))
        ! Remapping the pointer boundaries from Fortran default
        ! TODO(mrm): When supported use lower for c_f_pointer introduced in Fortran 2023
        call remap_fortran_pointer(tmat1, int([1, 1, mstart], kind=i32), int([mbsiz, ndim, mend], kind=i32))
        call remap_fortran_pointer(tmat2, int([1, 1, mstart], kind=i32), int([mbsiz, ndim, mend], kind=i32))

        ! We do not batch over frequencies as it will skyrocket the RAM consumption for the 
        ! tmat1 and tmat2
        OMP_OFFLOAD target data map(tofrom: wing1, wing2)
        do iom = iomstart, iomend
#if defined(FLANG_OPENMP_SLICE_MAP_BUG_WORKAROUND)
            OMP_OFFLOAD target has_device_addr(tmat1, tmat2) map(to: fnm)
#else
            OMP_OFFLOAD target has_device_addr(tmat1, tmat2) map(to: fnm(:,:,iom))
#endif
            !$omp teams distribute parallel do collapse(3) &
            !$omp default(none) private(ie2,ie1,imix) &
            !$omp shared(mstart,mend,ndim,mbsiz,tmat1,tmat2,fnm,minmmat,ik,iom)
            do ie2 = mstart, mend
                do ie1 = 1, ndim
                    do imix = 1, mbsiz
                        tmat1(imix,ie1,ie2) = fnm(ie1,ie2,iom)* &
                                                minmmat(imix,ie1,ie2)
                        tmat2(imix,ie1,ie2) = fnm(ie1,ie2,iom) * &
                                                conjg(minmmat(imix,ie1,ie2))
                    end do
                end do
            end do
            !$omp end teams distribute parallel do
            OMP_OFFLOAD end target

            ! Computing the wings here, note we do the three Cartesian
            ! directions in a single bunch, is almost an O2 operation
            ! i.e. memory bounded.
            ! We cannot batch it due to memory consumption
            call zgemm_gpu('n', 'n', mbsiz, 3, nmdim, coefw, tmat1_cptr, mbsiz, &
                        pm1_cptr, nmdim, zone, &
                        get_device_pointer(wing1(1,1,iom), my_device), mbsiz, device_world)
            call device_world%synchronize()
            call zgemm_gpu('n', 'n', mbsiz, 3, nmdim, coefw, tmat2_cptr, mbsiz, &
                        pm2_cptr, nmdim, zone, &
                        get_device_pointer(wing2(1,1,iom), my_device), mbsiz, device_world)
            call device_world%synchronize()

        end do
        OMP_OFFLOAD end target data

        ! Free memory
        nullify(pm1, pm2, tmat1, tmat2)

        call deallocate_device_memory(pm1_cptr, my_device)
        call deallocate_device_memory(pm2_cptr, my_device)
        call deallocate_device_memory(tmat1_cptr, my_device)
        call deallocate_device_memory(tmat2_cptr, my_device)

        call timesec(t1)
        time_dfwing = time_dfwing + (t1 - t0)
    end subroutine calcwings
end module mod_head_and_wings
