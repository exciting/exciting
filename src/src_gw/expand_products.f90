subroutine expand_products(ik,iq,nstart,nend,nsplit,mstart,mend,msplit,minm)

    use modmain, only : zzero
    use modgw,   only : mbsiz
    use precision, only: i32, dp
    use modgw,      only: kqset
    use iso_c_binding,         only: c_ptr, c_loc, c_f_pointer, c_sizeof, c_size_t
    use mod_pointer_remapping, only: remap_fortran_pointer
    use mod_device_offload,    only: device_world
    use device_linalg_common_interface, only: zgemm_gpu
    use m_memory_device,       only: allocate_device_memory, deallocate_device_memory, &
                                     bytes_double_complex, bytes_int, get_device_pointer
#include "offload.fpp"
    implicit none

    ! input variables
    integer(i32), intent(in) :: ik
    integer(i32), intent(in) :: iq
    integer(i32), intent(in) :: nstart, nend
    integer(i32), intent(in) :: nsplit         ! last index of valence states
    integer(i32), intent(in) :: mstart ,mend
    integer(i32), intent(in) :: msplit         ! last index of valence states
    complex(dp),  intent(inout) :: minm(mbsiz,nstart:nend,mstart:mend)

    ! local variables
    integer(i32) :: cstart, cend
    integer(i32) :: im, ie1, ie2
    integer(i32) :: iflag
    integer(i32) :: my_device
    integer(c_size_t) :: dim1, dim2, dim3
    type(c_ptr)  :: minm_valence_cptr
    complex(dp), pointer, contiguous :: minm_valence(:,:,:)
    type(c_ptr)  :: minm_core_cptr
    complex(dp), pointer, contiguous :: minm_core(:,:,:)

    dim1 = mbsiz
    my_device = device_world%get_device()

    OMP_OFFLOAD target
    ! IFX 2025.0.0 is not able to perform it efficiently except if a kernel is used
    ! Cray compiler generates a runtime error if using the kernel for the loop.
    ! GNU allows the latter so we go for it.
    ! For core states, we keep the Cray way for all cases, as the number
    ! of core states is always small in comparison with valence states
    ! and I want to keep compiler branching as small as possible.
    !
    ! TODO: Check future IFX compiler to see if this is improved.
    !
#if __INTEL_COMPILER
    !$omp teams distribute parallel do collapse(3) private(im,ie1,ie2)
    do ie2 = mstart, mend
      do ie1 = nstart, nend
        do im = 1, mbsiz
          minm(im,ie1,ie2) = zzero
        end do
      end do
    end do
    !$omp end teams distribute parallel do
#else
    minm(:,:,:) = zzero
#endif
    OMP_OFFLOAD end target

    if ((nsplit>0).and.(msplit<=0)) then

        !======================================================================
        ! When calculating the dielectric function:
        !
        ! n -> {1:nomax}{1:ncg}, so nsplit=nomax, when core states are treated
        ! m -> {numin:nstfv},       msplit is ignored
        !
        ! Important to specify msplit<=0
        ! to ignore the self energy related part (below)
        !
        !======================================================================
        if (nend<=nsplit) then

            ! calculate M^i_{nm}
            iflag = 1
            dim2  = nend - nstart + 1
            dim3  = mend   - mstart + 1
            call allocate_device_memory(minm_valence_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
            call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3])
            ! We are remapping the boundaries of a pointer with default Fortran lower bound
            ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
            !            this removes the need of remapping
            ! call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3], [1,nstart,mstart])
            ! The remapping is to prevent integer arithmetics in the loops
            call remap_fortran_pointer(minm_valence, int([1, nstart, mstart], kind=i32), int([mbsiz, nend, mend], kind=i32))

            call expand_products_block(ik,iq,nstart,nend,mstart,mend,minm_valence_cptr,iflag)

            OMP_OFFLOAD target has_device_addr(minm_valence)
#if __INTEL_COMPILER
            !$omp teams distribute parallel do collapse(3) private(im,ie1,ie2)
            do ie2 = mstart, mend
              do ie1 = nstart, nend
                do im = 1, mbsiz
                  minm(im,ie1,ie2) = minm_valence(im,ie1,ie2)
                end do
              end do
            end do
            !$omp end teams distribute parallel do
#else
            minm(:,nstart:nend,mstart:mend) = minm_valence(:,nstart:nend,mstart:mend)
#endif
            OMP_OFFLOAD end target

            nullify(minm_valence)
            call deallocate_device_memory(minm_valence_cptr, my_device)

        else if (nstart>nsplit) then

            ! calculate M^i_{cm}
            iflag = 2
            cstart = nstart-nsplit
            cend = nend-nsplit
            dim2 = cend - cstart + 1
            dim3 = mend - mstart + 1
            call allocate_device_memory(minm_core_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
            call c_f_pointer(minm_core_cptr, minm_core, [dim1,dim2,dim3])

            call expand_products_block(ik,iq,cstart,cend,mstart,mend,minm_core_cptr,iflag)

            OMP_OFFLOAD target has_device_addr(minm_core)
            minm(:,nstart:nend,mstart:mend) = minm_core(:,:,:)
            OMP_OFFLOAD end target

            nullify(minm_core)
            call deallocate_device_memory(minm_core_cptr, my_device)

        else

            ! calculate M^i_{nm}
            iflag = 1
            dim2  = nsplit - nstart + 1
            dim3  = mend   - mstart + 1
            call allocate_device_memory(minm_valence_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
            call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3])
            ! We are remapping the boundaries of a pointer with default Fortran lower bound
            ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
            !            this removes the need of remapping
            ! call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3], [1,nstart,mstart])
            ! The remapping is to prevent integer arithmetics in the loops
            call remap_fortran_pointer(minm_valence, int([1, nstart, mstart], kind=i32), int([mbsiz, nsplit, mend], kind=i32))

            call expand_products_block(ik,iq,nstart,nsplit,mstart,mend,minm_valence_cptr,iflag)

            OMP_OFFLOAD target has_device_addr(minm_valence)
#if __INTEL_COMPILER
            !$omp teams distribute parallel do collapse(3) private(im,ie1,ie2)
            do ie2 = mstart, mend
              do ie1 = nstart, nsplit
                do im = 1, mbsiz
                  minm(im,ie1,ie2) = minm_valence(im,ie1,ie2)
                end do
              end do
            end do
            !$omp end teams distribute parallel do
#else
            minm(:,nstart:nsplit,mstart:mend) = minm_valence(:,nstart:nsplit,mstart:mend)
#endif
            OMP_OFFLOAD end target

            nullify(minm_valence)
            call deallocate_device_memory(minm_valence_cptr, my_device)

            ! calculate M^i_{cm}
            iflag = 2
            cstart = 1
            cend = nend-nsplit
            dim2 = nend - nsplit - cstart + 1
            dim3  = mend   - mstart + 1
            call allocate_device_memory(minm_core_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
            call c_f_pointer(minm_core_cptr, minm_core, [dim1,dim2,dim3])

            call expand_products_block(ik,iq,cstart,cend,mstart,mend,minm_core_cptr,iflag)

            OMP_OFFLOAD target has_device_addr(minm_core)
            minm(:,nsplit+cstart:nend,mstart:mend) = minm_core(:,:,:)
            OMP_OFFLOAD end target

            nullify(minm_core)
            call deallocate_device_memory(minm_core_cptr, my_device)

        end if ! dielectric function part

    else if ((nsplit<=0).and.(msplit>0)) then

      !======================================================================
      ! When calculating the self energy:
      !
      ! n -> {ibgw:nbgw},      nsplit is ignored
      ! m -> {1:nstfv}{1:ncg}, so msplit=nstfv, when core states are treated
      !
      ! Important to specify nsplit<=0
      ! to ignore the dielectric function part (above)
      !
      !======================================================================
      if (mend<=msplit) then

        ! calculate M^i_{nm}
        iflag = 1
        dim2  = nend - nstart + 1
        dim3  = mend - mstart + 1
        call allocate_device_memory(minm_valence_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
        call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3])
        ! We are remapping the boundaries of a pointer with default Fortran lower bound
        ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
        !            this removes the need of remapping
        ! call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3], [1,nstart,mstart])
        ! The remapping is to prevent integer arithmetics in the loops
        call remap_fortran_pointer(minm_valence, int([1, nstart, mstart], kind=i32), int([mbsiz, nend, mend], kind=i32))

        call expand_products_block(ik,iq,nstart,nend,mstart,mend,minm_valence_cptr,iflag)

        OMP_OFFLOAD target has_device_addr(minm_valence)
#if __INTEL_COMPILER
        !$omp teams distribute parallel do collapse(3) private(im,ie1,ie2)
        do ie2 = mstart, mend
          do ie1 = nstart, nend
            do im = 1, mbsiz
              minm(im,ie1,ie2) = minm_valence(im,ie1,ie2)
            end do
          end do
        end do
        !$omp end teams distribute parallel do
#else
        minm(:,nstart:nend,mstart:mend) = minm_valence(:,nstart:nend,mstart:mend)
#endif
        OMP_OFFLOAD end target

        nullify(minm_valence)
        call deallocate_device_memory(minm_valence_cptr, my_device)

      else if (mstart>msplit) then

        ! calculate M^i_{nc}
        iflag = 3
        cstart = mstart-msplit
        cend = mend-msplit
        dim2 = nend - nstart + 1
        dim3 = cend - cstart + 1

        call allocate_device_memory(minm_core_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
        call c_f_pointer(minm_core_cptr, minm_core, [dim1,dim2,dim3])

        call expand_products_block(ik,iq,nstart,nend,cstart,cend,minm_core_cptr,iflag)

        OMP_OFFLOAD target has_device_addr(minm_core)
        minm(:,nstart:nend,mstart:mend) = minm_core(:,:,:)
        OMP_OFFLOAD end target

        nullify(minm_core)
        call deallocate_device_memory(minm_core_cptr, my_device)

      else
        ! calculate M^i_{nm}
        iflag = 1
        dim2  = nend   - nstart + 1
        dim3  = msplit - mstart + 1
        call allocate_device_memory(minm_valence_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
        call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3])
        ! We are remapping the boundaries of a pointer with default Fortran lower bound 
        ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
        !            this removes the need of remapping
        ! call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3], [1,nstart,mstart])
        ! The remapping is to prevent integer arithmetics in the loops
        call remap_fortran_pointer(minm_valence, int([1, nstart, mstart], kind=i32), int([mbsiz, nend, msplit], kind=i32))

        call expand_products_block(ik,iq,nstart,nend,mstart,msplit,minm_valence_cptr,iflag)

        OMP_OFFLOAD target has_device_addr(minm_valence)
#if __INTEL_COMPILER
        !$omp teams distribute parallel do collapse(3) private(im,ie1,ie2)
        do ie2 = mstart, msplit
          do ie1 = nstart, nend
            do im = 1, mbsiz
              minm(im,ie1,ie2) = minm_valence(im,ie1,ie2)
            end do
          end do
        end do
        !$omp end teams distribute parallel do
#else
        minm(:,nstart:nend,mstart:msplit) = minm_valence(:,nstart:nend,mstart:msplit)
#endif
        OMP_OFFLOAD end target

        nullify(minm_valence)
        call deallocate_device_memory(minm_valence_cptr, my_device)

        ! calculate M^i_{nc}
        iflag = 3
        cstart = 1
        cend = mend-msplit
        dim2 = nend - nstart + 1
        dim3 = mend - msplit - cstart + 1
        call allocate_device_memory(minm_core_cptr, dim1*dim2*dim3*bytes_double_complex, my_device)
        call c_f_pointer(minm_core_cptr, minm_core, [dim1,dim2,dim3])

        call expand_products_block(ik,iq,nstart,nend,cstart,cend,minm_core_cptr,iflag)

        OMP_OFFLOAD target has_device_addr(minm_core)
        minm(:,nstart:nend,msplit+cstart:mend) = minm_core(:,:,:)
        OMP_OFFLOAD end target

        nullify(minm_core)
        call deallocate_device_memory(minm_core_cptr, my_device)

      end if ! self energy part

    else

      write(*,*) 'emergency stop (expand_product): wrong usage of the subroutine!'
      write(*,*) 'currently implemented options are:'
      write(*,*) "nsplit>0 and msplit<=0 used when calculating the dielectric function"
      write(*,*) "nsplit<=0 and msplit>0 used when calculating the self energy"
      stop

    end if

    return

contains

    !=============================================================================
    !
    !=============================================================================
    subroutine expand_products_block(ik,iq,nstart,nend,mstart,mend,minm_cptr,iflag)
        use modmain,               only: zone, zzero
        use modgw,                 only: locmatsiz, matsiz, mbsiz, fgw
        use mod_coulomb_potential, only: barc
        use device_linalg_common_interface, only: zgemm_gpu
        use mod_mpi_gw
        implicit none
        ! input variables
        integer(i32), intent(in) :: ik    ! the index of the first k-point
        integer(i32), intent(in) :: iq    ! the index of the q-point
        integer(i32), intent(in) :: nstart, nend  ! range of n states
        integer(i32), intent(in) :: mstart, mend  ! range of m states
        type(c_ptr),  value      :: minm_cptr ! The C pointer to the minm we are working on
        integer(i32), intent(in) :: iflag ! 1-calcminm; 2-calcmicm; 3-calcminc
        ! local variables
        integer(i32) :: nmdim
        complex(dp), allocatable :: minm_(:,:,:)

        nmdim = (nend - nstart + 1) * (mend - mstart + 1)

        select case(iflag)

          case(1)
            ! Valence contribution
            ! calculate v^{1/2}*M^i_{nm}
            allocate(minm_(matsiz,nstart:nend,mstart:mend))
            OMP_OFFLOAD target data map(alloc: minm_)
            call calcminm2(ik,iq,nstart,nend,mstart,mend,minm_)
            call zgemm_gpu('c','n', &
                       mbsiz,nmdim,matsiz, &
                       zone, &
                       get_device_pointer(barc,my_device),matsiz, &
                       get_device_pointer(minm_,my_device),matsiz, &
                       zzero,minm_cptr,mbsiz,device_world)
            call device_world%synchronize()
            OMP_OFFLOAD end target data
            deallocate(minm_)

          case(2)
            ! Core contribution
            ! calculate v^{1/2}*M^i_{cm}
            allocate(minm_(locmatsiz,nstart:nend,mstart:mend))
            call calcmicm(ik,iq,nstart,nend,mstart,mend,minm_)
            OMP_OFFLOAD target data map(alloc: minm_)
            call zgemm_gpu('c','n', &
                       mbsiz,nmdim,locmatsiz, &
                       zone, &
                       get_device_pointer(barc,my_device),matsiz, &
                       get_device_pointer(minm_,my_device),locmatsiz, &
                       zzero,minm_cptr,mbsiz,device_world)
            call device_world%synchronize()
            OMP_OFFLOAD end target data
            deallocate(minm_)

          case(3)
            ! Core contribution
            ! calculate v^{1/2}*M^i_{nc}
            allocate(minm_(locmatsiz,nstart:nend,mstart:mend))
            call calcminc(ik,iq,nstart,nend,mstart,mend,minm_)
            OMP_OFFLOAD target data map(alloc: minm_)
            call zgemm_gpu('c','n', &
                       mbsiz,nmdim,locmatsiz, &
                       zone, &
                       get_device_pointer(barc,my_device),matsiz, &
                       get_device_pointer(minm_,my_device),locmatsiz, &
                       zzero,minm_cptr,mbsiz,device_world)
            call device_world%synchronize()
            OMP_OFFLOAD end target data
            deallocate(minm_)

          case default
            write(*,*) 'ERROR(expand_products:expand_products_block):'
            write(*,*) 'Unknown iflag=', iflag

        end select

      return
    end subroutine

end subroutine



