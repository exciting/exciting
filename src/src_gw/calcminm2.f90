!!-----------------------------------------------------------------------------------
!! This routine calculates the expansion coefficients, i.e. M^i_nm(k,q) (see Eqs. 50 and 53 of 10.1016/j.cpc.2012.09.018),
!! for n and m being valence and conduction states.
!!-----------------------------------------------------------------------------------
subroutine calcminm2(ik,iq,nstart,nend,mstart,mend,minm)
#if _CRAYFTN
! Instruct the Cray compiler to use aggressive optimization for this file
!DIR$ OPTIMIZE(-haggress)
#endif
  use modinput
  use modmain,                only : nspecies, natmtot, natoms, idxas, idxlm, idxlo, &
                                    intgv, apword, nlorb, lorbl, &
                                    nlomax, apwordmax, nmatmax, nlotot
  use constants,              only : zzero, zone, pi, twopi, zi
  use modgw,                  only : kqset, Gkqset, Gqbarc, Gqset, Gset, fdebug, time_minm
  use mod_bands,              only : eveck, eveckp, eveckalm, eveckpalm
  use mod_product_basis,      only : nmix, bigl, bradketa, bradketlo, mpwipw, &
                                     matsiz, locmatsiz, mbindex
  use mod_misc_gw,            only : vi, atposl
  use mod_gaunt_coefficients, only : epsangint, getgauntcoef 
  use iso_c_binding,          only : c_ptr, c_loc, c_f_pointer, c_sizeof, c_size_t
  use mod_device_offload,     only : device_world
  use device_linalg_common_interface, only : zgemm_batched_gpu
  use m_memory_device,        only : allocate_device_memory, deallocate_device_memory, &
                                     bytes_double_complex, bytes_int, get_device_pointer
  use mod_pointer_remapping,  only : remap_fortran_pointer 
  use precision,              only : i32, dp
#include "offload.fpp"

  implicit none
  integer(i32), intent(in) :: ik
  integer(i32), intent(in) :: iq
  integer(i32), intent(in) :: nstart, nend
  integer(i32), intent(in) :: mstart, mend
  complex(dp), intent(out):: minm(matsiz,nstart:nend,mstart:mend)

  integer(i32) :: jk
  integer(i32) :: bl, bm
  integer(i32) :: i, ia, is, ias
  integer(i32) :: igk1, igk2
  integer(i32) :: idxlo1, idxlo2
  integer(i32) :: io1, io2, ilo1, ilo2
  integer(i32) :: imix, irm
  integer(i32) :: ie1, ie2
  integer(i32) :: l1, m1, l2, m2, l1m1, l2m2
  integer(i32) :: l2min, l2max
  integer(i32) :: ig
  integer(i32), dimension(3):: ikv, ig0 ! Indexes of G_1+G'-G
  integer(i32) :: ngk1, ngk2
  integer(i32) :: ndim, mdim
  integer(i32) :: l_max_apw
  integer(i32) :: imix_start, imix_end

  real(dp) :: x, arg, angint
  real(dp) :: qvec(3)
  real(dp) :: tstart, tend, tmt

  real(dp) :: t1, t2

  complex(dp) :: phs, bk, sqvi
  complex(dp), allocatable :: lok(:,:), lokp(:,:)
  type(c_ptr) :: tmat_cptr, veckp_cptr
  complex(dp), contiguous, pointer :: tmat(:,:,:), veckp(:,:)
  complex(dp), allocatable :: phase(:)

  ! Variables for the offload of the interstitial part
  type(c_ptr) :: igqk12_cptr
  integer(i32), contiguous, pointer :: igqk12(:,:)
  type(c_ptr) :: temp_int_cptr, temp2_int_cptr, mnn_cptr, eveck_int_cptr, eveckp_int_cptr
  complex(dp), contiguous, pointer :: temp_int(:,:,:), temp2_int(:,:,:), mnn(:,:,:), eveck_int(:,:),  eveckp_int(:,:)

  integer(i32) :: nblocks_ngq, ngq_block_counter, igq_start, igq_end, igq_block, ngq_max_block_size, igq, ngq

  integer(i32) :: my_device

  call timesec(tstart)

  ! Get some offload data
  my_device = device_world%get_device()

  ! GPU-aware compilation: minm is assumed to be in the device
  ! IFX 2025.0.0 is not able to perform it efficiently except if a kernel is used
  ! Cray compiler generates a runtime error if using the kernel for the loop
  ! GNU allows the latter so we go for it.
  ! For core states, we keep the Cray way for all cases, as the number
  ! of core states is always small in comparison with valence states
  ! and I want to keep compiler branching as small as possible.
  !
  ! TODO: Check future IFX compiler to see if this is improved.
  !
  OMP_OFFLOAD target
#if __INTEL_COMPILER
  !$omp teams distribute parallel do collapse(3) default(none) &
  !$omp shared(matsiz,nstart,nend,mstart,mend,minm) private(imix,ie1,ie2)
  do ie2 = mstart, mend
    do ie1 = nstart, nend
      do imix = 1, matsiz
        minm(imix,ie1,ie2) = zzero
      end do
    end do
  end do
  !$omp end teams distribute parallel do
#else
  minm(:,:,:) = zzero
#endif
  OMP_OFFLOAD end target

  ! index ranges
  ndim = nend-nstart+1
  mdim = mend-mstart+1
  l_max_apw = input%groundstate%lmaxapw

  jk = kqset%kqid(ik,iq)
  do i = 1, 3
    qvec(i) = kqset%vql(i,iq)
    x = kqset%vkl(i,ik)-kqset%vkl(i,jk)-qvec(i)
    ig0(i) = nint(x)
  end do

  !======================
  ! MT region
  !======================

  ! Copy in transposed way the local orbitals
  ! cofficients for improved memory access.

  allocate(lok(nstart:nend,nlotot))
  do ie1 = nstart, nend
    lok(ie1,1:)=eveck(Gkqset%ngk(1,ik)+1:Gkqset%ngk(1,ik)+nlotot,ie1)
  end do

  allocate(lokp(mstart:mend,nlotot))
  do ie1 = mstart,mend
    lokp(ie1,1:) = eveckp(Gkqset%ngk(1,jk)+1:Gkqset%ngk(1,jk)+nlotot,ie1)
  end do

  ! Computing the phase for the Mixed Basis - MuffinTin, i.e. (\e^{-iqr_\alpha}\)
  ! where (\\alpha\) is the atomic index. This is done here
  ! as exp is an slow function the call of which should be minimized,
  ! especially in the device.
  allocate(phase(natmtot))
  do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia,is)
      arg = twopi * dot_product(atposl(1:3,ia,is), qvec(1:3))
      phase(ias) = exp(- cmplx(0.0_dp, arg, dp))
    end do
  end do

  ! Create a temporary array for the MT part of the expansion coefficient with an improved
  ! memory pattern, so that the instensive loop can be deployed as an efficient kernel
  ! to the device.
  call allocate_device_memory(tmat_cptr, ndim * mdim * locmatsiz * bytes_double_complex, my_device)
  call c_f_pointer(tmat_cptr, tmat, int([ndim, mdim, locmatsiz], kind=c_size_t))
  ! We are remapping the boundaries of a pointer with default Fortran lower bound 
  ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
  !            this removes the need of remapping
  ! call c_f_pointer(tmat_cptr, tmat, int([ndim, mdim, locmatsiz], kind=c_size_t), &
  !                  int([nstart, mstart, 1], kind=c_size_t))
  ! The remapping is to prevent integer arithmetics in the loops
  call remap_fortran_pointer(tmat, int([nstart, mstart, 1],kind=i32), int([nend, mend, locmatsiz],kind=i32))

  ! Allocate a temporary array to store the vectors of the outer products
  ! Note that in devices (GPUs), that "big" allocatable private arrays are not behaving very well, as it tries to allocate in the
  ! shared memory of the CU, which is quite small. Therefore, the second rank is to allow for it to be "shared", as a unique "copy"
  ! exist for each mixed basis, which prevents race conditions.
  call allocate_device_memory(veckp_cptr, mdim * locmatsiz * bytes_double_complex, my_device)
  call c_f_pointer(veckp_cptr, veckp, int([mdim, locmatsiz], kind=c_size_t))
  ! We are remapping the boundaries of a pointer with default Fortran lower bound
  ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
  !            this removes the need of remapping
  ! call c_f_pointer(veckp_cptr, veckp, int([mdim, locmatsiz], kind=c_size_t), &
  !                  int([mstart,1], kind=c_size_t))
  ! The remapping is to prevent integer arithmetics in the loops
  call remap_fortran_pointer(veckp, int([mstart, 1],kind=i32), int([mend, locmatsiz],kind=i32))

  OMP_OFFLOAD target enter data map(always, to: lok, lokp, phase)

  !> Computing now the actual MT contribution
  !> \[
  !>   M_{n m}^i(\mathbf{k}, \mathbf{q}) = e^{-i \mathbf{q} \cdot \mathbf{r}_\alpha}
  !>    \sum_{l_1 m_1} \sum_{l_2 m_2}\left[g_{L M ; l_2 m_2}^{l_1 m_1}\right]^*
  !>    \sum_{\zeta_1 \zeta_2} \mathcal{A}_{n \mathbf{k}, \alpha \zeta_1 l_1 m_1}
  !>    \mathcal{A}_{m \mathbf{k}-\mathbf{q}, \alpha \zeta_2 l_2 m_2}^*
  !>    \left\langle N L \mid \zeta_1 l_1, \zeta_2 l_2\right\rangle_\alpha
  !> \]
  !
  !  Note(mrm): When modifying this loop, avoid private allocatables as they are 
  !             not properly handled by the GPUs. Also, keep innermost loops over 
  !             memory contiguous slices.

  ! For GPU aware compilation

  ! First partition our workspace so that teams process groups of mixed functions
  ! Teams are assigned tasks accessing non-contiguous memory, as they run in different
  ! excution units sharing different fast memory.
  ! The threads of each team will take care of teams internal loops, in particular those involving bands, which
  ! consume most of the time (indeed the outer products consume 80% of the loop time in the CPU). These operate on
  ! contiguous data, as they share fast access memory. Those are invoked with OMP_OFFLOAD parallel do

  ! For CPU-only compilation the mixed basis in the MT are scattered across OpenMP threads

  OMP_OFFLOAD target has_device_addr(tmat, veckp)
  OMP_OFFLOAD teams distribute &
  OMP_NO_OFFLOAD parallel do schedule(dynamic) &
  !$omp default(none) private(imix,ie2,ie1,is,ia,irm,bl,bm,ias,arg,phs,l1,m1,l1m1,l2min,l2max) &
  !$omp private(l2,m2,l2m2,angint,io1,io2,ilo2,idxlo1,ilo1,idxlo2,bk) &
  !$omp shared(locmatsiz,mstart,mend,nstart,nend,mbindex,idxas,l_max_apw,idxlm,apword,bradketa,eveckpalm) &
  !$omp shared(nlorb,lorbl,idxlo,lokp,eveckalm,bradketlo,lok,phase,tmat,veckp,ndim)
  do imix = 1, locmatsiz

    ! Getting species id, atomic id for given is, the radial function id, the L of the MB, and the "m" of the MB
    is  = mbindex(imix,1)
    ia  = mbindex(imix,2)
    irm = mbindex(imix,3)
    bl  = mbindex(imix,4)
    bm  = mbindex(imix,5)

    ! Atomic idx and phase
    ias = idxas(ia,is)
    phs = phase(ias)

    ! Set tmat for the given mixed basis to zero. 
    ! OMP_OFFLOAD parallel do: In device divide the job among the threads within a execution unit.
    ! In the rest of the imix loop these sections mean the
    ! same parallelization level.
    OMP_OFFLOAD parallel do schedule(static,1) collapse(2)
    do ie2 = mstart, mend
      do ie1 = nstart, nend
        tmat(ie1,ie2,imix) = zzero
      end do
    end do
    OMP_OFFLOAD end parallel do

    ! Sum over l1m1 and l2m2
    ! Notice the index of the summation
    do l1 = 0, l_max_apw
      do m1 = -l1, l1
        l1m1 = idxlm(l1,m1) !! sum over l1m1

        l2min = abs(bl-l1)
        l2max = min(bl+l1,l_max_apw)
        do l2 = l2min, l2max
          m2 = -bm+m1
          if (abs(m2) <= l2) then
            l2m2 = idxlm(l2,m2) !! sum over l2m2

            ! Angular integral. (Gaunt coefficient)
            angint = getgauntcoef(l2,bl,l1,m2,bm)
            if (abs(angint) < epsangint) cycle

            do io1 = 1, apword(l1,is) !! sum over \zeta_1
              !======
              ! APW-APW
              !======
              OMP_OFFLOAD parallel do schedule(static,1)
              do ie2 = mstart, mend
                veckp(ie2, imix) = zzero
              end do
              OMP_OFFLOAD end parallel do

              do io2 = 1, apword(l2,is) !! sum over \zeta_2
                bk = angint*bradketa(2,irm,l1,io1,l2,io2,ias)
                OMP_OFFLOAD parallel do schedule(static,1)
                do ie2 = mstart, mend
                  veckp(ie2, imix)=veckp(ie2, imix)+bk*eveckpalm(ie2,io2,l2m2,ias)
                end do
                OMP_OFFLOAD end parallel do
              end do
              !======
              ! APW-LO
              !======
              do ilo2 = 1, nlorb(is)  !! sum over local orbital (LO_2)
                if (lorbl(ilo2,is)==l2) then
                  idxlo2 = idxlo(l2m2,ilo2,ias)
                  bk = angint*bradketa(3,irm,l1,io1,ilo2,1,ias)
                  OMP_OFFLOAD parallel do schedule(static,1)
                  do ie2 = mstart, mend
                    veckp(ie2,imix)=veckp(ie2,imix)+bk*lokp(ie2,idxlo2)
                  end do
                  OMP_OFFLOAD end parallel do
                end if
              end do ! ilo2

              ! Outer product (Computationally intensive part of the loop)
              OMP_OFFLOAD parallel do collapse(2) schedule(static,1)
              do ie2 = mstart, mend
                do ie1 = nstart, nend
                  tmat(ie1,ie2,imix) = tmat(ie1,ie2,imix) + eveckalm(ie1,io1,l1m1,ias) * veckp(ie2,imix)
                end do ! ie1
              end do ! ie2
              OMP_OFFLOAD end parallel do

            end do ! io1


            do ilo1 = 1, nlorb(is)   !! sum over local orbital (LO_1)
              if (lorbl(ilo1,is)==l1) then
                idxlo1 = idxlo(l1m1,ilo1,ias)
                !======
                ! LO-A
                !======
                OMP_OFFLOAD parallel do schedule(static,1)
                do ie2 = mstart, mend
                  veckp(ie2, imix) = zzero
                end do
                OMP_OFFLOAD end parallel do
                do io2 = 1, apword(l2,is) !! sum over \zeta_2
                  bk = angint*bradketlo(2,irm,ilo1,l2,io2,ias)
                  OMP_OFFLOAD parallel do schedule(static,1)
                  do ie2 = mstart, mend
                    veckp(ie2,imix)=veckp(ie2,imix)+bk*eveckpalm(ie2,io2,l2m2,ias)
                  end do
                  OMP_OFFLOAD end parallel do
                end do ! io2
                !======
                ! LO-LO
                !======
                do ilo2 = 1, nlorb(is) !! sum over local orbital (LO_2)
                  if (lorbl(ilo2,is)==l2) then
                    idxlo2 = idxlo(l2m2,ilo2,ias)
                    bk = angint*bradketlo(3,irm,ilo1,ilo2,1,ias)
                    OMP_OFFLOAD parallel do schedule(static,1)
                    do ie2 = mstart, mend
                      veckp(ie2,imix)=veckp(ie2,imix)+bk*lokp(ie2,idxlo2)
                    end do
                    OMP_OFFLOAD end parallel do
                  end if
                end do ! ilo2

                ! Outer product (Computationally intensive part of the loop)
                OMP_OFFLOAD parallel do collapse(2) schedule(static,1)
                do ie2 = mstart, mend
                  do ie1 = nstart, nend
                    tmat(ie1,ie2,imix) = tmat(ie1,ie2,imix) + lok(ie1,idxlo1) * veckp(ie2,imix)
                  end do ! ie1
                end do ! ie2
                OMP_OFFLOAD end parallel do

              end if
            end do ! ilo1

          end if ! m2
        end do ! l2

      end do ! m1
    end do ! l1

    OMP_OFFLOAD parallel do schedule(static,1) collapse(2)
    do ie2 = mstart, mend
      do ie1 = nstart, nend
        tmat(ie1,ie2,imix) = phs * tmat(ie1,ie2,imix)
      end do ! ie1
    end do ! ie2
    OMP_OFFLOAD end parallel do

  end do ! imix
  OMP_OFFLOAD end teams distribute
  OMP_OFFLOAD end target
  OMP_NO_OFFLOAD end parallel do


  ! Now copy the the matrix in the proper order
  OMP_OFFLOAD target has_device_addr(tmat)
  !$omp teams distribute parallel do collapse(3) &
  !$omp shared(locmatsiz,nstart,nend,mstart,mend,tmat,minm) &
  !$omp private(imix,ie1,ie2)
  do imix = 1, locmatsiz
    do ie2 = mstart, mend
      do ie1 = nstart, nend
        minm(imix,ie1,ie2) = tmat(ie1,ie2,imix)
      end do
    end do
  end do
  !$omp end teams distribute parallel do
  OMP_OFFLOAD end target

  ! Free memory
  nullify(tmat, veckp)
  call deallocate_device_memory(veckp_cptr, my_device)
  call deallocate_device_memory(tmat_cptr, my_device)

  OMP_OFFLOAD target exit data map(delete: lokp, lok, phase)
  deallocate(lok, lokp, phase)

  !======================
  ! Interstitial region (Eq. 53)
  !======================
  ! Loop over the mixed basis functions:
  sqvi = cmplx(sqrt(vi),0.0_dp,kind=dp)
  ngk1 = Gkqset%ngk(1,ik)
  ngk2 = Gkqset%ngk(1,jk)

  call allocate_device_memory(igqk12_cptr, ngk1 * ngk2 * bytes_int, my_device)
  call c_f_pointer(igqk12_cptr, igqk12, [int(ngk1,kind=c_size_t), int(ngk2,kind=c_size_t)])

  OMP_OFFLOAD target map(to: ig0) has_device_addr(igqk12)
  !$omp teams distribute parallel do collapse(2) &
  !$omp default(NONE) shared(ik,jk,iq,ngk2,ngk1,Gset,Gkqset,ig0) &
  !$omp shared(igqk12,Gqbarc) &
  !$omp private(igk2,igk1,ikv,ig)
  do igk2 = 1, ngk2 ! loop over G
    do igk1 = 1, ngk1 ! loop over G'
      ikv(1:3) = Gset%ivg(1:3,Gkqset%igkig(igk1,1,ik)) - &
                 Gset%ivg(1:3,Gkqset%igkig(igk2,1,jk)) + ig0(1:3)
      if( all( ikv >= Gset%intgv(:, 1) ) .and. all( ikv <= Gset%intgv(:, 2) ) ) then
        ig = Gset%ivgig(ikv(1),ikv(2),ikv(3))
        igqk12(igk1,igk2) = Gqbarc%igigk(ig,1,iq)
      else
        igqk12(igk1,igk2) = -1
      end if
    end do ! igk2
  end do ! igk1
  !$omp end teams distribute parallel do
  OMP_OFFLOAD end target

  ngq_max_block_size = merge(Gqset%ngk(1,iq), min(input%gw%GBatchCount, Gqset%ngk(1,iq)), input%gw%GBatchCount <= 0)
  nblocks_ngq        = ceiling(real(Gqset%ngk(1,iq), kind = dp) /  ngq_max_block_size)

  call allocate_device_memory(eveck_int_cptr, ngk1 * ndim * bytes_double_complex, my_device)
  call allocate_device_memory(eveckp_int_cptr, ngk2 * mdim * bytes_double_complex, my_device)
  call allocate_device_memory(temp_int_cptr, ngk2 * ngk1 * ngq_max_block_size * bytes_double_complex, my_device)
  call allocate_device_memory(temp2_int_cptr, ngk2 * ndim * ngq_max_block_size * bytes_double_complex, my_device)
  call allocate_device_memory(mnn_cptr, mdim * ndim * ngq_max_block_size * bytes_double_complex, my_device)

  call c_f_pointer(eveck_int_cptr, eveck_int, int([ngk1,ndim],c_size_t))
  call c_f_pointer(eveckp_int_cptr, eveckp_int, int([ngk2,mdim],c_size_t))

  OMP_OFFLOAD target has_device_addr(eveck_int, eveckp_int)
  eveck_int(:,:)  = eveck(1:ngk1, nstart:nend)
  eveckp_int(:,:) = eveckp(1:ngk2, mstart:mend)
  OMP_OFFLOAD end target

  ! We iterate over batches of plane waves, because for 
  ! systems with very large number of planes waves operating 
  ! over the full set can obliterate the RAM.
  do ngq_block_counter = 1, nblocks_ngq

    igq_start = (ngq_block_counter - 1) * ngq_max_block_size + 1
    igq_end   = min(Gqset%ngk(1,iq), igq_start + ngq_max_block_size - 1)
    ngq       = igq_end - igq_start + 1

    call c_f_pointer(temp_int_cptr, temp_int, int([ngk2,ngk1,ngq],c_size_t))
    call c_f_pointer(temp2_int_cptr, temp2_int, int([ngk2,ndim,ngq],c_size_t))
    call c_f_pointer(mnn_cptr, mnn, int([mdim, ndim, ngq],c_size_t))

    ! We are remapping the boundaries of a pointer with default Fortran lower bound
    ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
    !            this removes the need of remapping
    ! call c_f_pointer(temp_int_cptr, temp_int, int([ngk2,ngk1,ngq],c_size_t), int([1, 1, igq_start]))
    ! call c_f_pointer(mnn_cptr, mnn, int([mdim, ndim, ngq],c_size_t), int([1, 1, igq_start]))
    ! The remapping is to prevent integer arithmetics in the loops
    call remap_fortran_pointer(temp_int, int([1, 1, igq_start],kind=i32), int([ngk2, ngk1, igq_end],kind=i32))
    call remap_fortran_pointer(mnn, int([mstart,nstart,igq_start],kind=i32), int([mend, nend, igq_end],kind=i32))

    OMP_OFFLOAD target has_device_addr(igqk12,temp_int)
    !$omp teams distribute parallel do collapse(3) &
    !$omp default(NONE) private(igq,igk1,igk2) &
    !$omp shared(ngk1,ngk2,mpwipw,temp_int,igqk12,igq_start,igq_end)
    do igq = igq_start, igq_end
      do igk1 = 1, ngk1 ! loop over G+k
        do igk2 = 1, ngk2 ! loop over G'+k'
          if (igqk12(igk1,igk2) > 0) then
            temp_int(igk2,igk1,igq) = mpwipw(igq,igqk12(igk1,igk2))
          else
            temp_int(igk2,igk1,igq) = zzero
          end if
        end do ! igk2
      end do ! igk1
    end do
    !$omp end teams distribute parallel do
    OMP_OFFLOAD end target

    ! Matrix-matrix product either in the device or in the host depending on the compilation
    ! NVIDIA and AMD: calls magma zgemm
    ! Intel: calls MKL device offloaded version
    ! CPU: call the zgemm of the base library
    ! 
    ! Note(mrm): For efficiency purposes the zgemms are batched
    !
    call zgemm_batched_gpu('n', 'n', ngk2, ndim, ngk1, zone, temp_int_cptr, ngk2, &
                  eveck_int_cptr, ngk1, zzero, temp2_int_cptr, ngk2, ngq, device_world, strideb=0_i32)
    call device_world%synchronize()
    call zgemm_batched_gpu('t', 'n', mdim, ndim, ngk2, sqvi, eveckp_int_cptr, ngk2, &
                  temp2_int_cptr, ngk2, zzero, mnn_cptr, mdim, ngq, device_world, stridea=0_i32)
    call device_world%synchronize()

    OMP_OFFLOAD target has_device_addr(mnn)
    !$omp teams distribute parallel do collapse(3) &
    !$omp default(NONE) private(igq,ie2,ie1) shared(mstart,mend,nstart,nend,locmatsiz,igq_start,igq_end,minm,mnn)
    do igq = igq_start, igq_end
      do ie2 = mstart, mend
        do ie1 = nstart, nend
          minm(locmatsiz+igq,ie1,ie2) = mnn(ie2,ie1,igq)
        end do ! ie2
      end do ! ie1
    end do ! igq
    !$omp end teams distribute parallel do
    OMP_OFFLOAD end target

    nullify(temp_int, temp2_int, mnn)

  end do ! ngq_block_counter

  ! Free C-allocated memory
  nullify(eveck_int, eveckp_int, igqk12)
  call deallocate_device_memory(eveck_int_cptr, my_device)
  call deallocate_device_memory(eveckp_int_cptr, my_device)
  call deallocate_device_memory(temp_int_cptr, my_device)
  call deallocate_device_memory(temp2_int_cptr, my_device)
  call deallocate_device_memory(mnn_cptr, my_device)
  call deallocate_device_memory(igqk12_cptr, my_device)

  call timesec(tend)
  time_minm = time_minm + tend - tstart

end subroutine calcminm2

