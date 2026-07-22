
subroutine calcmwm(ikp, jk, nstart, nend, mstart, mend, minm, offdiagonal)
    use constants, only: zone, zzero, pi
    use modgw,   only: vi, kqset, Gamma, singc1, singc2, mbsiz, &
                       epsilon, epsh, epsw1, epsw2, mwm, freq, freq_selfc
    use precision, only: dp, i32
    use mod_pointer_remapping, only: remap_fortran_pointer 
    use mod_device_offload,    only: device_world
    use device_linalg_common_interface, only: zgemm_batched_gpu, zhemm_gpu, zgemm_gpu, zhemv_gpu
    use m_memory_device,       only: allocate_device_memory, deallocate_device_memory, &
                                   bytes_double_complex, bytes_int, get_device_pointer
    use iso_c_binding, only: c_ptr, c_size_t, c_f_pointer
    use mod_atoms, only: idxas
    use mod_core_states, only: corind
    use mod_bands, only: evalfv, nstse
    use mod_corestate, only: evalcr
    use mod_eigenvalue_occupancy, only: efermi
    use mod_offdiagonal_selfenergy, only: add_q_omega_contrib_to_offdiagonal_selfenergy_correl_at_ik

#include "offload.fpp"

    implicit none

    integer(i32),  intent(in) :: ikp
    integer(i32), intent(in)  :: jk
    integer(i32), intent(in) :: nstart, nend
    integer(i32), intent(in) :: mstart, mend
    complex(dp), intent(in)  :: minm(mbsiz, nstart:nend, mstart:mend)
    logical, intent(in)      :: offdiagonal

    ! local variables
    integer(i32) :: iom
    integer(i32) :: ndim, mdim, nmdim, nomeg
    integer(i32) :: ie1, ie2, iemin, iemax, ieg
    integer(i32) :: my_device
    real(dp)     :: vi4pi, wkq
    complex(dp)  :: coefs1, coefs2
    type(c_ptr)  :: wm_cptr
    complex(dp), contiguous, pointer :: wm(:,:,:)

    my_device = device_world%get_device()

    wkq    = 1.0_dp/real(kqset%nkpt,kind=dp)
    ndim   = nend - nstart + 1
    mdim   = mend - mstart + 1
    nmdim  = mdim*ndim
    nomeg  = size( epsilon, 3 )

    if (Gamma) then
      vi4pi  = 4.0_dp*pi*vi
      coefs1 = singc1*sqrt(vi4pi)
      coefs2 = singc2*vi4pi
      iemin  = max(nstart,mstart)
      iemax  = min(nend,mend)
    end if 

    !-------------------------------------------------
    ! calculate \sum_{ij} M^i_{nm}* W^c_{ij} M^j_{nm}
    !-------------------------------------------------
    call allocate_device_memory(wm_cptr, int(mbsiz*nmdim,kind=c_size_t)*bytes_double_complex, my_device)
    call c_f_pointer(wm_cptr, wm, int([mbsiz,ndim,mdim],kind=c_size_t))
    ! We are remapping the boundaries of a pointer with default Fortran lower bound 
    ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
    !            this removes the need of remapping
    ! call c_f_pointer(wm_cptr, wm, int([mbsiz,ndim,mdim],kind=c_size_t), int([1,nstart,mstart],kind=c_size_t))
    ! The remapping is to prevent integer arithmetics in the loops
    call remap_fortran_pointer(wm, int([1, nstart, mstart], kind=i32), int([mbsiz, nend, mend], kind=i32))

    do iom = 1, nomeg

      ! Compute W^c_{ij} M^j_{nm} for a given frequency
      call  zhemm_gpu('l','u',mbsiz,nmdim,zone,get_device_pointer(epsilon(1,1,iom), my_device), mbsiz, &
                      get_device_pointer(minm(1,nstart,mstart), my_device), mbsiz, &
                      zzero, wm_cptr, mbsiz, device_world)
      call device_world%synchronize()

      ! Performing the batch of dot_products in parallel. 
      OMP_OFFLOAD target has_device_addr(wm)
      !$omp teams distribute parallel do collapse(2) default(none) private(ie2, ie1)&
      !$omp shared(mstart, mend, nstart, nend, mwm, wkq, minm, wm, iom)
      do ie2 = mstart, mend
        do ie1 = nstart, nend
          mwm(ie1,ie2,iom) = wkq*dot_product(minm(:,ie1,ie2),wm(:,ie1,ie2))
        end do
      end do
      !$omp end teams distribute parallel do
      OMP_OFFLOAD end target

#if !defined(FLANG_OPENMP_SLICE_MAP_BUG_WORKAROUND)
      OMP_OFFLOAD target update from(mwm(nstart:nend,mstart:mend,iom))
#endif

      ! Compute the contribution of the MWM product to the offdiagonal
      ! correlation self-energy.
      if (offdiagonal) then
        call add_q_omega_contrib_to_offdiagonal_selfenergy_correl_at_ik( & 
               nstart, nend, mstart, mend, nstse, evalfv, &
               evalcr, corind, idxas, efermi, &
               minm, wm, wkq, freq, freq_selfc, &
               nomeg, iom, ikp, jk)
      end if

    end do ! iom
 
#if defined(FLANG_OPENMP_SLICE_MAP_BUG_WORKAROUND)
      OMP_OFFLOAD target update from(mwm)
#endif

    nullify(wm)
    call deallocate_device_memory(wm_cptr, my_device)

    ! Add the contribution of the head and the wings
    ! Notice that in device that is not the most efficient, as the
    ! memory access pattern is not contiguous. 
    if (Gamma) then
      !$omp parallel do collapse(2) default(none) &
      !$omp private(ieg, iom) &
      !$omp shared(mwm, coefs1, coefs2, minm, epsh, epsw2, epsw1, mbsiz, iemin, iemax, nomeg)
      do iom = 1, nomeg
        do ieg = iemin, iemax
          mwm(ieg,ieg,iom) = mwm(ieg,ieg,iom) + &
                            coefs2*epsh(1,1,iom) + &
                            coefs1*(dot_product(conjg(minm(:,ieg,ieg)),epsw2(:,1,iom)) + &
                                    dot_product(minm(:,ieg,ieg),epsw1(:,1,iom)))
        end do
      end do
      !$omp end parallel do 

#if defined(FLANG_OPENMP_SLICE_MAP_BUG_WORKAROUND)
      OMP_OFFLOAD target update to(mwm)
#endif

    end if

    return
end subroutine

