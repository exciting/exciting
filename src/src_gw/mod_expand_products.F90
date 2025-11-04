!> Module providing expansion coefficients \( M^i_{nm} \)
!> that are used to expand products of two Kohn-Sham wavefunctions
!> in either the mixed basis or in the v-diagonal basis.
!> n and m can each be either valence or core states.
!> All four possible combinations of valence and core
!> are organized in blocks.
!>```
!>         -- m -->
!>   +--------+--------+
!>   |        |        |
!>   |   VV   |   VC   |
!> | |        |        |
!> n +--------+--------+
!> | |        |        |
!> v |   CV   |   CC   |
!>   |        |        |
!>   +--------+--------+
!>```
!> Note: core-core contribution is not yet implemented
module mod_expand_products
  use modmain,                        only: zzero, zone
  use modgw,                          only: mbsiz, matsiz, locmatsiz
  use precision,                      only: i32, dp
  use iso_c_binding,                  only: c_ptr, c_f_pointer, c_size_t
  use mod_pointer_remapping,          only: remap_fortran_pointer
  use mod_device_offload,             only: device_world
  use device_linalg_common_interface, only: zgemm_gpu
  use m_memory_device,                only: allocate_device_memory,   &
                                            deallocate_device_memory, &
                                            bytes_double_complex,     &
                                            get_device_pointer
  use mod_coulomb_potential,          only: barc
  use modmpi,                         only: terminate_if_false
  use xstring,                        only: newline
  use asserts,                        only: assert
  implicit none

  !> Flag indicating computation type: valence-valence 
  integer(i32), parameter :: FLAG_VV=1
  !> Flag indicating computation type: core-valence
  integer(i32), parameter :: FLAG_CV=2
  !> Flag indicating computation type: valence-core
  integer(i32), parameter :: FLAG_VC=3
  !> Flag indicating computation type: core-core
  integer(i32), parameter :: FLAG_CC=4

  private
  public :: expand_products_generic, split_interval
#include "offload.fpp"

contains

  !> Main routine for expanding wavefunction products into their valence/core
  !> blocks. This routine initializes the output array and delegates the actual
  !> computations to specialized helper routines.
  subroutine expand_products_generic( ik, iq,                   &
                                      n_val_start, n_val_end,   &
                                      n_core_start, n_core_end, &
                                      m_val_start, m_val_end,   &
                                      m_core_start, m_core_end, &
                                      minm,                     &
                                      apply_barc )
    !> k-point index
    integer(i32), intent(in)                  :: ik
    !> q-point index
    integer(i32), intent(in)                  :: iq
    !> first valence-n index
    integer(i32), intent(in)                  :: n_val_start
    !> last  valence-n index
    integer(i32), intent(in)                  :: n_val_end
    !> first core-n index
    integer(i32), intent(in)                  :: n_core_start
    !> last  core-n index
    integer(i32), intent(in)                  :: n_core_end
    !> first valence-m index
    integer(i32), intent(in)                  :: m_val_start
    !> last  valence-m index
    integer(i32), intent(in)                  :: m_val_end
    !> first core-m index
    integer(i32), intent(in)                  :: m_core_start
    !> last  core-m index
    integer(i32), intent(in)                  :: m_core_end
    !> Array of expansion coefficients in either
    !> mixed basis or v-diagonal (Coulomb) basis
    !> The array's 2nd dim runs from
    !> n_val_start…n_val_end + n_core_start…n_core_end,
    !> and 3rd dim likewise covers
    !> m_val_start…m_val_end + m_core_start…m_core_end.
    complex(dp),  intent(inout), allocatable  :: minm(:,:,:)
    !> If true, store \( v^{1/2} M^i_{nm} \) in the v-diagonal (Coulomb-eigen)
    !> basis; otherwise store \( M^i_{nm} \) in the mixed basis.
    logical,      intent(in)                  :: apply_barc

    ! local variables
    integer(i32)                              :: nstart,nend,mstart,mend
    integer(i32)                              :: ie1,ie2,im
    integer(i32)                              :: dim_n, dim_n_val, dim_n_core
    integer(i32)                              :: dim_m, dim_m_val, dim_m_core
    integer(i32)                              :: basis_size
    integer(i32)                              :: n_core_offset, m_core_offset

    call assert( allocated(minm), 'minm must be allocated before call')
    call assert( lbound(minm,1) == 1, &
         'First dimension of minm must have lower bound 1.' )
    
    ! Determine the basis size. It differs between mbsiz and matsiz depending on the 
    ! application or not of the bare Coulomb potential
    basis_size = size(minm,1)

    if (apply_barc) then
       call assert(allocated(barc), 'barc must be allocated if apply_barc=.true.')
       call assert( basis_size == mbsiz, &
            'First dimension of minm must equal mbsiz when v-diagonal basis (barc) is applied.' )
    else
       call assert( basis_size == matsiz, &
            'First dimension of minm must equal matsiz when barc is not applied.' )
    end if

    ! Determine minm bounds
    nstart = lbound(minm,2); nend   = ubound(minm,2)
    mstart = lbound(minm,3); mend   = ubound(minm,3)

    dim_n      = nend - nstart + 1
    dim_n_val  = max(n_val_end - n_val_start + 1, 0_i32)
    dim_n_core = max(n_core_end - n_core_start + 1, 0_i32)
    dim_m      = mend - mstart + 1
    dim_m_val  = max(m_val_end - m_val_start + 1, 0_i32)
    dim_m_core = max(m_core_end - m_core_start + 1, 0_i32)

    call assert( dim_n_val + dim_n_core == dim_n, &
         'Dimension mismatch for n in M^i_{nm}' )
    call assert( dim_m_val + dim_m_core == dim_m, &
         'Dimension mismatch for m in M^i_{nm}' )
    call assert( (dim_n_val==0) .or. (lbound(minm,2) <= n_val_start .and. n_val_end <= ubound(minm,2)), &
         'Valence n-subrange must lie within minm second-dimension bounds' )
    call assert( (dim_m_val==0) .or. (lbound(minm,3) <= m_val_start .and. m_val_end <= ubound(minm,3)), &
         'Valence m-subrange must lie within minm third-dimension bounds' )
    
    ! Zero full array (Intel-optimized if available)
    OMP_OFFLOAD target
#if __INTEL_COMPILER
    !$omp teams distribute parallel do collapse(3)
    do ie2 = mstart, mend; do ie1 = nstart, nend; do im = 1,basis_size
       minm(im,ie1,ie2) = zzero
    end do; end do; end do
    !$omp end teams distribute parallel do
#else
    minm(:,:,:) = zzero
#endif
    OMP_OFFLOAD end target

    ! Compute offset for placing core blocks inside the global minm layout
    n_core_offset = block_offset( n_core_start, n_core_end, &
                                  n_val_start,  n_val_end,  &
                                  lbound(minm,2), ubound(minm,2) )
    m_core_offset = block_offset( m_core_start, m_core_end, &
                                  m_val_start,  m_val_end,  &
                                  lbound(minm,3), ubound(minm,3) )
    
    ! Dispatch blocks if non-empty
    if(n_val_start<=n_val_end .and. m_val_start<=m_val_end) then
       call expand_products_block_wrapper(ik,iq,n_val_start,n_val_end,m_val_start,m_val_end, &
            FLAG_VV,minm, 0_i32, 0_i32, apply_barc)
    end if
    if(n_core_start<=n_core_end .and. m_val_start<=m_val_end) then
       call expand_products_block_wrapper(ik,iq,n_core_start,n_core_end,m_val_start,m_val_end, &
            FLAG_CV,minm, n_core_offset, 0_i32, apply_barc)
    end if
    if(n_val_start<=n_val_end .and. m_core_start<=m_core_end) then
       call expand_products_block_wrapper(ik,iq,n_val_start,n_val_end,m_core_start,m_core_end, &
            FLAG_VC,minm, 0_i32, m_core_offset, apply_barc)
    end if
    if(n_core_start<=n_core_end .and. m_core_start<=m_core_end) then
       call expand_products_block_wrapper(ik,iq,n_core_start,n_core_end,m_core_start,m_core_end, &
            FLAG_CC,minm, n_core_offset, m_core_offset, apply_barc)
    end if
  end subroutine expand_products_generic

  !> Helper: allocates device buffer, maps slice, calls block, copies back
  subroutine expand_products_block_wrapper( ik, iq, n1, n2, m1, m2, flag, minm, n_offset, m_offset, apply_barc )
    !> k-point index
    integer(i32), intent(in)                 :: ik
    !> q-point index
    integer(i32), intent(in)                 :: iq
    !> first n index
    integer(i32), intent(in)                 :: n1
    !> last n index
    integer(i32), intent(in)                 :: n2
    !> first m index
    integer(i32), intent(in)                 :: m1
    !> last m index
    integer(i32), intent(in)                 :: m2
    !> Computation flag for valence/core combination (VV, CV, VC, CC)
    integer(i32), intent(in)                 :: flag
    !> Array of expansion coefficients same as in generic (with full coverage in the 2nd/3rd dims)
    complex(dp),  intent(inout), allocatable :: minm(:,:,:)
    !> shift to apply to the 2nd index of minm
    integer(i32), intent(in)                 :: n_offset
    !> shift to apply to the 3rd index of minm
    integer(i32), intent(in)                 :: m_offset
    !> If true, store \( v^{1/2} M^i_{nm} \) in the v-diagonal (Coulomb-eigen)
    !> basis; otherwise store \( M^i_{nm} \) in the mixed basis.
    logical,      intent(in)                 :: apply_barc

    integer(i32)  :: d1, d2, d3
    integer(i32)  :: dev, ie1, ie2, im
    type(c_ptr)   :: buf_cptr
    complex(dp), pointer, contiguous :: buf(:,:,:)
    integer(c_size_t) :: num_of_elements
    integer(i32)  :: nloc, nrows
    
    dev = device_world%get_device()

    d1 = size(minm,1)
    d2 = n2-n1+1
    d3 = m2-m1+1
    num_of_elements = int(d1, kind=c_size_t) * int(d2, kind=c_size_t) * int(d3, kind=c_size_t)


    call allocate_device_memory(buf_cptr, num_of_elements * bytes_double_complex, dev)
    call c_f_pointer(buf_cptr,buf,[d1,d2,d3])
    ! We are remapping the boundaries of a pointer with default Fortran lower bound
    ! TODO(mrm): In Fortran2023 the lower is added to CALL C_F_POINTER(cptr, fptr [, shape, lower])
    !            this removes the need of remapping
    ! call c_f_pointer(minm_valence_cptr, minm_valence, [dim1,dim2,dim3], [1,nstart,mstart])
    ! The remapping is to prevent integer arithmetics in the loops
    call remap_fortran_pointer(buf, int([1,n1,m1],i32), int([d1,n2,m2],i32))

    call expand_products_block(ik, iq, n1, n2, m1, m2, buf_cptr, flag, apply_barc)

    nloc = merge(matsiz, locmatsiz, flag==FLAG_VV)
    nrows = merge(mbsiz, nloc, apply_barc)  ! if barc: mbsiz, else nloc
    
    OMP_OFFLOAD target has_device_addr(buf)
#if __INTEL_COMPILER
    !$omp teams distribute parallel do collapse(3)
    do ie2 = m1,m2; do ie1 = n1,n2; do im = 1,nrows
       minm(im,n_offset+ie1,m_offset+ie2) = buf(im,ie1,ie2)
    end do; end do; end do
    !$omp end teams distribute parallel do
#else
    minm(1:nrows,n_offset+n1:n_offset+n2,m_offset+m1:m_offset+m2) = buf(1:nrows,n1:n2,m1:m2)
#endif
    OMP_OFFLOAD end target

    nullify(buf)
    call deallocate_device_memory(buf_cptr, dev)
  end subroutine expand_products_block_wrapper

  !> Performs the low-level computation for a specific block type
  !> (valence-valence, core-valence, valence-core, core-core). 
  subroutine expand_products_block( ik, iq, nstart, nend, mstart, mend, buf_cptr, flag, apply_barc )
    !> k-point index
    integer(i32), intent(in)  :: ik
    !> q-point index
    integer(i32), intent(in)  :: iq
    !> first n index
    integer(i32), intent(in)  :: nstart
    !> last n index
    integer(i32), intent(in)  :: nend
    !> first m index
    integer(i32), intent(in)  :: mstart
    !> last m index
    integer(i32), intent(in)  :: mend
    !> Device pointer to the output buffer
    type(c_ptr), value        :: buf_cptr
    !> Computation flag for valence/core combination (VV, CV, VC, CC)
    integer(i32), intent(in)  :: flag
    !> If true, store \( v^{1/2} M^i_{nm} \) in the v-diagonal (Coulomb-eigen)
    !> basis; otherwise store \( M^i_{nm} \) in the mixed basis.
    logical,    intent(in)    :: apply_barc

    ! local variables
    complex(dp), allocatable  :: tmp(:,:,:)
    integer(i32)              :: nmdim, nloc, ie1, ie2, im
    integer(i32)              :: dev
    integer, parameter :: msglen = 200
    character(len=msglen) :: error_message
    complex(dp), pointer, contiguous :: buf(:,:,:)

    dev   = device_world%get_device()
    nmdim = (nend-nstart+1)*(mend-mstart+1)

    nloc = merge(matsiz, locmatsiz, flag==FLAG_VV)

    allocate(tmp(nloc,nstart:nend,mstart:mend))
    OMP_OFFLOAD target data map(alloc:tmp)

    select case(flag)
    case(FLAG_VV)
       ! Valence-Valence
       ! For GPU compilation tmp is filled in the GPU
       ! so no transfer is needed
       call calcminm2(ik,iq,nstart,nend,mstart,mend,tmp)
    case(FLAG_CV)
       ! Core-Valence
       call calcmicm(ik,iq,nstart,nend,mstart,mend,tmp)
       ! Uploading to the GPU
       OMP_OFFLOAD target update to(tmp)
    case(FLAG_VC)
       ! Valence-Core
       call calcminc(ik,iq,nstart,nend,mstart,mend,tmp)
       ! Upload to the GPU
       OMP_OFFLOAD target update to(tmp)
    case(FLAG_CC)
       write(error_message, '(A,A,A)') &
            'ERROR(expand_products:expand_products_block):', newline, &
            'Core-Core contribution is not yet implemented.'
       call terminate_if_false(.false., trim(error_message))
       ! Core-Core
       ! call calcmicc(ik,iq,nstart,nend,mstart,mend,tmp)
    case default
       write(error_message, '(A,A,A,I0)') &
            'ERROR(expand_products:expand_products_block):', newline, &
            'Unknown flag=', flag
       call terminate_if_false(.false., trim(error_message))
    end select

    if (apply_barc) then
       ! Multiply by barc:  v^{1/2} * tmp -> buf
       call zgemm_gpu('c','n', mbsiz, nmdim, nloc, &
            zone, get_device_pointer(barc,dev), matsiz, &
            get_device_pointer(tmp,dev),        nloc, &
            zzero, buf_cptr, mbsiz, device_world)
    else
       ! TODO(mrm): We need a better solution, this is wasting memory
       ! and resources doing a copy
       call c_f_pointer(buf_cptr, buf, shape(tmp))
       call remap_fortran_pointer(buf, lbound(tmp), ubound(tmp))
       
       OMP_OFFLOAD target has_device_addr(buf)
#if __INTEL_COMPILER
       !$omp teams distribute parallel do collapse(3)
       do ie2 = mstart,mend; do ie1 = nstart,nend; do im = 1,nloc
          buf(im,ie1,ie2) = tmp(im,ie1,ie2)
       end do; end do; end do
       !$omp end teams distribute parallel do
#else
       buf(:,nstart:nend,mstart:mend) = tmp(:,nstart:nend,mstart:mend)
#endif
       OMP_OFFLOAD end target
       nullify(buf)
    end if
    call device_world%synchronize()
    OMP_OFFLOAD end target data

    deallocate(tmp)
  end subroutine expand_products_block


  !> Splits a given interval into valence and core regions based on a split point.
  !> On exit, any region with start>end should be treated as empty.
  subroutine split_interval(     &
       mstart, mend, msplit,     &
       m_val_start, m_val_end,   &
       m_core_start, m_core_end)

    !-- Inputs
    !> lower bound of original interval
    integer(i32), intent(in)  :: mstart
    !> upper bound of original interval
    integer(i32), intent(in)  :: mend
    !> core/valence splitting point
    integer(i32), intent(in)  :: msplit        

    !-- Outputs
    !> start of valence-region
    integer(i32), intent(out) :: m_val_start
    !> end of valence-region
    integer(i32), intent(out) :: m_val_end
    !> start of core-region (core indices after msplit start with 1 again)
    integer(i32), intent(out) :: m_core_start
    !> end of core-region (core indices after msplit start with 1 again)
    integer(i32), intent(out) :: m_core_end   

    !---- valence-region: [mstart, min(mend, msplit)] ----------------------
    m_val_start = mstart
    m_val_end   = min(mend, msplit)
    ! If m_val_end < m_val_start, this region is empty (start>end).

    !---- core-region: [mstart-msplit, mend-msplit], clamped for emptiness --
    m_core_start = max(mstart - msplit, 1)
    m_core_end   = max(mend   - msplit, 0)
    ! If m_core_end < m_core_start, this region is empty (start>end).

  end subroutine split_interval
  
  !> Compute offset for placing the core block inside the global minm layout.
  !>
  !> This routine provides robust handling of all combinations of valence
  !> and core blocks with arbitrary Fortran lower bounds. The following
  !> rules are applied:
  !>
  !>   * If both valence and core blocks are present:
  !>     Pack core immediately after valence, i.e.
  !>     shifted_core_start = valence_end + 1
  !>     (No gaps between blocks.)
  !>
  !>   * If the valence block is empty (core-only case):
  !>     Preserve original core positions relative to the global bounds
  !>     such that:
  !>     shifted_core_start = global_lbound
  !>     (No artificial index shift.)
  !>
  !>   * If the core block is empty:
  !>     Return offset = 0  (unused).
  !>
  !> On exit, the returned offset may be added to each core-index:
  !> ie1_shifted = core_index + offset
  integer(i32) function block_offset( core_start, core_end,  &
                              valence_start,  valence_end,           &
                              global_lbound, global_ubound )     &
                              result(offset)
    !-- Inputs
    !> first core index in original indexing
    integer(i32), intent(in) :: core_start
    !> last core index in original indexing
    integer(i32), intent(in) :: core_end
    !> first valence index in original indexing
    integer(i32), intent(in) :: valence_start
    !> last valence index in original indexing
    integer(i32), intent(in) :: valence_end
    !> global lower bound (e.g. lbound(minm,2))
    integer(i32), intent(in) :: global_lbound
    !> global upper bound (e.g. ubound(minm,2))
    integer(i32), intent(in) :: global_ubound

    !-- Local flags
    logical :: has_core, has_valence

    has_core = (core_start <= core_end)
    has_valence  = (valence_start  <= valence_end)

    if (has_valence) then
       ! Pack core directly behind valence:  (valence_end+1) → first core
       offset = (valence_end + 1_i32) - core_start
    else
       ! Core-only: preserve original positioning relative to global bounds
       !   i.e. shifted_core_start == global_lbound
       offset = global_lbound - core_start
    end if

    if (has_core) then
       !-- Safety check: shifted core block must lie within global bounds
       call assert( (offset + core_start) >= global_lbound .and. &
            (offset + core_end)   <= global_ubound,              &
            'block_offset: shifted core block out of global bounds' )
    else
       ! Empty core block → offset unused
       offset = 0_i32
    end if

  end function block_offset

  
end module mod_expand_products
