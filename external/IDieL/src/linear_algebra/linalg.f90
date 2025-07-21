! Copyright (C) 2020-2024 GreenX library
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied. See the License for the specific language governing
! permissions and limitations under the License.

!> @file
!> This file contain subroutines containing linear algebra using MAGMA (GPU) or LAPACK (CPU) like library
!> the routines are specifically written down case by case so the GPU sending and retrieving 
!> is minimized

!> This module contains the subroutines that build up vectors and blocks using linear algebra
module idiel_linalg

    use idiel_constants, only: i32, i32, aip, zzero, zone, fourpi
    use iso_c_binding


! Set some macros so the proper linear algebra are called
#ifdef USE_SINGLE_PRECISION
#define _getrf_gpu         cgetrf_gpu
#define _getri_gpu         cgetri_gpu
#define _get_getri_nb_gpu  get_cgetri_nb_gpu
#define _gemm_gpu          cgemm_gpu
#define _sizecomplex       bytes_single_complex
#else
#define _getrf_gpu         zgetrf_gpu
#define _getri_gpu         zgetri_gpu
#define _get_getri_nb_gpu  get_zgetri_nb_gpu
#define _gemm_gpu          zgemm_gpu
#define _sizecomplex       bytes_double_complex
#endif   

    use device_linalg_common_interface, only: _getrf_gpu, _getri_gpu, _get_getri_nb_gpu, _gemm_gpu
    use m_memory_device, only: get_device_pointer, _sizecomplex, allocate_device_memory, &
                               deallocate_device_memory

    use m_device_world_t, only: device_world_t
#include "offload.fpp"

    implicit none

contains

    !> Inverts the A matrix and save the result to inverse_A
    !> @param[in] A - the complex matrix to invert
    !> @param[in] inverse_A - the inverse of A
    !> @param[in] world - the gpu_world_t handler
    subroutine inverse_complex_LU(A, inverse_A, world)

        complex(aip), contiguous, intent(in)            :: A(:,:)
        complex(aip), target, allocatable, intent(out)  :: inverse_A(:,:)
        type(device_world_t), intent(inout)             :: world

        ! LAPACK
        integer(i32) :: info, lwork, nb, n
        integer(i32), allocatable :: ipiv(:)
        type(C_ptr) :: dA_ptr, dwork_ptr

        ! Some constants
        n = size(A,1)
        ! Allocate
        allocate(inverse_A, source = A)

        nb = _get_getri_nb_gpu(n)
        lwork = nb * size(A,1)
        allocate(ipiv(size(A,1)))

        ! For device compilation
        OMP_OFFLOAD target data map(tofrom: inverse_A)
        ! Allocate in device work array
        call allocate_device_memory(dwork_ptr, _sizecomplex * lwork, world%get_device())

        ! Obtain device pointers
        dA_ptr    = get_device_pointer(inverse_A, world%get_device())
        
        call _getrf_gpu(n, n, dA_ptr, n, ipiv, info, world)
        if (info /= 0) error stop "inverse_complex_LU: error calling Xgetrf gpu"
        call _getri_gpu(n, dA_ptr, n, ipiv, dwork_ptr, lwork, info, world)
        if (info /= 0) error stop "inverse_complex_LU: error calling Xgetri gpu"

        call deallocate_device_memory(dwork_ptr, world%get_device())
        OMP_OFFLOAD end target data

        deallocate(ipiv)

    end subroutine inverse_complex_LU

    !> It computes the auxiliary S and T vectors and the macroscopic dielectric matrix
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] S       - the auxiliary vector related to lower wing
    !> @param[in] wingU    - the upper wing
    !> @param[out] T       - the auxiliary vector related to upper wing
    !> @param[inout] L     - the macroscopic dielectric ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager   
    subroutine compute_auxiliary_and_macroscopic_3d_general(Binv, wingL, S, wingU, T, L, world)

        complex(aip), target, intent(inout)                        :: Binv(:,:)
        complex(aip), target, intent(inout)                        :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)             :: S(:,:)
        complex(aip), target, intent(inout)                        :: wingU(:,:)
        complex(aip), target, allocatable, intent(out)             :: T(:,:)
        complex(aip), target, intent(inout)                        :: L(:,:)
        type(device_world_t), intent(inout)                        :: world

        integer(i32) :: nb
        type(C_ptr) :: L_dptr, S_dptr, T_dptr, Binv_dptr, wingL_dptr, wingU_dptr

        nb = size(Binv,1)
        allocate(S(nb,3))
        allocate(T(3,nb))
        
        OMP_OFFLOAD target enter data map(to: Binv, wingL, wingU, L)
        OMP_OFFLOAD target enter data map(alloc: S, T)

        ! Get device pointers
        L_dptr     = get_device_pointer(L, world%get_device())
        S_dptr     = get_device_pointer(S, world%get_device())
        T_dptr     = get_device_pointer(T, world%get_device())
        Binv_dptr  = get_device_pointer(Binv, world%get_device())
        wingL_dptr = get_device_pointer(wingL, world%get_device())
        wingU_dptr = get_device_pointer(wingU, world%get_device())

        call _gemm_gpu('n', 'n', nb, 3, nb, zone, Binv_dptr, &
                    nb, wingL_dptr, nb, zzero, S_dptr, nb, world)
        call _gemm_gpu('t', 'n', 3, nb, nb, zone, wingU_dptr, &
                    nb, Binv_dptr, nb, zzero, T_dptr, 3, world)
        call _gemm_gpu('t', 'n', 3, 3, nb, -zone, wingU_dptr, &
                    nb, S_dptr, nb, zone, L_dptr, 3, world)
        call world%synchronize()
        
        ! Update L, S, T
        OMP_OFFLOAD target update from(L, S, T)

        ! Delete from device
        OMP_OFFLOAD target exit data map(delete: Binv, wingL, wingU, L, S, T)

    end subroutine compute_auxiliary_and_macroscopic_3d_general


    !> It computes the auxiliary S auxiliary vector and the macroscopic dielectric matrix
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] S       - the auxiliary vector related to lower wing
    !> @param[inout] L     - the macroscopic dielectric ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager
    subroutine compute_auxiliary_and_macroscopic_3d_hermitian(Binv, wingL, S, L, world)

        complex(aip), target, intent(inout)                        :: Binv(:,:)
        complex(aip), target, intent(inout)                        :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)             :: S(:,:)
        complex(aip), target, intent(inout)                        :: L(:,:)
        type(device_world_t), intent(inout)                        :: world

        integer(i32) :: nb

        type(C_ptr) :: L_dptr, S_dptr, T_dptr, Binv_dptr, wingL_dptr

        nb = size(Binv,1)
        allocate(S(nb,3))

        OMP_OFFLOAD target enter data map(to: Binv, wingL, L)
        OMP_OFFLOAD target enter data map(alloc: S)

        ! Get device pointers
        L_dptr     = get_device_pointer(L, world%get_device())
        S_dptr     = get_device_pointer(S, world%get_device())
        Binv_dptr  = get_device_pointer(Binv, world%get_device())
        wingL_dptr = get_device_pointer(wingL, world%get_device())

        call _gemm_gpu('n', 'n', nb, 3, nb, zone, Binv_dptr, &
                        nb, wingL_dptr, nb, zzero, S_dptr, nb, world)
        call _gemm_gpu('c', 'n', 3, 3, nb, -zone, wingL_dptr, &
                        nb, S_dptr, nb, zone, L_dptr, 3, world)
        call world%synchronize()

        ! Deassociate from device but keep them there
        ! Update L, S
        OMP_OFFLOAD target update from(L, S)

        ! Delete from device
        OMP_OFFLOAD target exit data map(delete: Binv, wingL, L, S)

    end subroutine compute_auxiliary_and_macroscopic_3d_hermitian

    !> It computes 1 / q^T \cdot L \cdot q 
    !> @param[in] L      - macroscopic dielectric matrix (if GPU it must be loaded)
    !> @param[in] q      - the object q (if GPU it must be loaded)
    !> @param[in] world  - the gpu_world_t handler 
    !> @returns invqLq (if GPU it will also provide the GPU object)
    subroutine compute_inverse_head(L, q, world, invqLq)

        complex(aip),    target, allocatable, intent(inout) :: L(:,:)
        complex(aip),    target, allocatable, intent(in)    :: q(:,:)
        type(device_world_t), intent(inout)                 :: world
        complex(aip),    target, allocatable, intent(out)   :: invqLq(:)

        complex(aip), pointer, contiguous :: Lq(:,:)
        integer(i32) :: nr 
        integer(i32) :: i

        type(C_ptr) :: q_dptr, Lq_dptr, L_dptr 
        integer     :: device

        nr = size(q,2)

        call allocate_device_memory(Lq_dptr, _sizecomplex * product(shape(q)), world%get_device())
        call c_f_pointer(Lq_dptr, Lq, shape(q))

        OMP_OFFLOAD target data map(to: L)

        ! Get device pointers
        L_dptr   = get_device_pointer(L, world%get_device())
        q_dptr   = get_device_pointer(q, world%get_device())

        call _gemm_gpu('n', 'n', 3, nr, 3, zone, L_dptr, &
                         3, q_dptr, 3, zzero, Lq_dptr, 3, world)
        call world%synchronize()

        if (allocated(invqLq)) deallocate(invqLq)
        allocate(invqLq(size(q,2)))

        OMP_OFFLOAD target map(tofrom: invqLq) has_device_addr(Lq)
        !$omp teams distribute parallel do private(i)
        do i = 1, size(q,2)
                invqLq(i) = 1.0_aip / dot_product(q(:,i),Lq(:,i))
        end do
        !$omp end teams distribute parallel do
        OMP_OFFLOAD end target

        nullify(Lq)
        call deallocate_device_memory(Lq_dptr, world%get_device())

        OMP_OFFLOAD end target data

        deallocate(L)

    end subroutine compute_inverse_head


    !> It computes Sq term of  \frac{ \mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}}
    !> @param[in] q - the linear algebra object containing q
    !> @param[in] S - S vector
    !> @param[in] world - the gpu_world_t handler 
    !> @param[out] qS   - (in GPU if compiled for)
    subroutine compute_inverse_wingL(q, S, world, qS)
        complex(aip), target, allocatable,  intent(inout)   :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)   :: S(:,:)
        type(device_world_t),       intent(inout)           :: world
        complex(aip), target, allocatable,  intent(inout)   :: qS(:,:)

        integer(i32) :: nr, nb 

        type(C_ptr) :: qS_dptr, S_dptr, q_dptr

        nr = size(q,2)
        nb = size(S,1)
        allocate(qS(nr,nb))

        OMP_OFFLOAD target enter data map(to: S)
        OMP_OFFLOAD target enter data map(alloc: qS)

        ! Get device pointers
        q_dptr   = get_device_pointer(q, world%get_device())
        qS_dptr  = get_device_pointer(qS, world%get_device())
        S_dptr   = get_device_pointer(S, world%get_device())

        call _gemm_gpu('t', 't', nr, nb, 3, zone, q_dptr, &
                        3, S_dptr, nb, zzero, qS_dptr, nr, world)
        call world%synchronize()

        OMP_OFFLOAD target update from(qS)
        OMP_OFFLOAD target exit data map(delete: S, qS)

        deallocate(S)

    end subroutine compute_inverse_wingL

    !> It computes Tq term of \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}}
    !> @param[in] ref_q - the linear algebra object containing q
    !> @param[in] ref_S - the linear algebra object containing T
    !> @param[in] world - the gpu_world_t handler 
    !> @param[out] qT   - (in GPU if compiled for)
    subroutine compute_inverse_wingU(q, T, world, qT)
        complex(aip), target, allocatable,  intent(inout)  :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)  :: T(:,:)
        type(device_world_t),               intent(inout)  :: world
        complex(aip), target, allocatable,  intent(out)    :: qT(:,:)

        integer(i32) :: nr, nb 

        type(C_ptr) :: qT_dptr, T_dptr, q_dptr

        nr = size(q,2)
        nb = size(T,2)
        allocate(qT(nr,nb))

        OMP_OFFLOAD target enter data map(to: T)
        OMP_OFFLOAD target enter data map(alloc: qT)

        ! Get device pointers
        q_dptr   = get_device_pointer(q, world%get_device())
        qT_dptr  = get_device_pointer(qT, world%get_device())
        T_dptr   = get_device_pointer(T, world%get_device())

        call _gemm_gpu('t', 'n', nr, nb, 3, zone, q_dptr, &
                        3, T_dptr, 3, zzero, qT_dptr, nr, world)
        
        call world%synchronize()
        OMP_OFFLOAD target update from(qT)
        OMP_OFFLOAD target exit data map(delete: qT)

        deallocate(T)

    end subroutine compute_inverse_wingU

    !> It computes the auxiliary ag, bg vectors and the A tensor
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] ag      - the auxiliary vector related to lower wing
    !> @param[in] wingU    - the upper wing
    !> @param[out] bg      - the auxiliary vector related to upper wing
    !> @param[inout] A     - the A ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager   
    subroutine compute_auxiliary_and_A_2d_general(Binv, wingL, ag, wingU, bg, A, world)
        complex(aip), target, intent(inout)              :: Binv(:,:)
        complex(aip), target, intent(inout)              :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)   :: ag(:,:)
        complex(aip), target, intent(inout)              :: wingU(:,:)
        complex(aip), target, allocatable, intent(out)   :: bg(:,:)
        complex(aip), target, intent(inout)              :: A(:,:)
        type(device_world_t), intent(inout)                 :: world

        integer :: nb, i

        type(C_ptr) :: A_dptr, ag_dptr, bg_dptr, Binv_dptr, wingL_dptr, wingU_dptr

        nb = size(Binv,1)
        allocate(ag(nb,3))

        ! Set arrays to the ones required by the formalism
        do i = 1, 3
            A(i,i) = zone - A(i,i)
        end do

        wingL = -cmplx(1.0_aip/sqrt(fourpi), 0.0_aip, aip) * wingL
        A     = -cmplx(1.0_aip/fourpi, 0.0_aip, aip) * A
        wingU = -cmplx(1.0_aip/sqrt(fourpi), 0.0_aip, aip) * wingU
        allocate(bg(3,nb))

        OMP_OFFLOAD target enter data map(to: Binv, wingL, wingU, A)
        OMP_OFFLOAD target enter data map(alloc: ag, bg)

        ! Get device pointers
        A_dptr     = get_device_pointer(A    , world%get_device())
        ag_dptr    = get_device_pointer(ag   , world%get_device())
        bg_dptr    = get_device_pointer(bg   , world%get_device())
        Binv_dptr  = get_device_pointer(Binv , world%get_device())
        wingL_dptr = get_device_pointer(wingL, world%get_device())
        wingU_dptr = get_device_pointer(wingU, world%get_device())

        call _gemm_gpu('n', 'n', nb, 3, nb, -zone, Binv_dptr, &
                    nb, wingL_dptr, nb, zzero, ag_dptr, nb, world)
        call _gemm_gpu('t', 'n', 3, nb, nb, -zone, wingU_dptr, &
                    nb, Binv_dptr, nb, zzero, bg_dptr, 3, world)
        call _gemm_gpu('t', 'n', 3, 3, nb, -zone, wingU_dptr, &
                    nb, ag_dptr, nb, zone, A_dptr, 3, world)
        call world%synchronize()

        OMP_OFFLOAD target update from(A, ag, bg)
        OMP_OFFLOAD target exit data map(delete: Binv, wingL, wingU, A, ag, bg)

    end subroutine compute_auxiliary_and_A_2d_general

    !> It computes the auxiliary ag vector and the A tensor
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] ag      - the auxiliary vector related to lower wing
    !> @param[inout] A     - the A ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager
    subroutine compute_auxiliary_and_A_2d_hermitian(Binv, wingL, ag, A, world)
        complex(aip), target, intent(inout)                   :: Binv(:,:)
        complex(aip), target, intent(inout)                   :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)        :: ag(:,:)
        complex(aip), target, intent(inout)                   :: A(:,:)
        type(device_world_t), intent(inout)                      :: world

        integer :: nb, i

        type(C_ptr) :: A_dptr, ag_dptr, bg_dptr, Binv_dptr, wingL_dptr

        nb = size(Binv,1)
        allocate(ag(nb,3))

        ! Set arrays to the ones required by the formalism
        do i = 1, 3
            A(i,i) = zone - A(i,i)
        end do

        wingL = -cmplx(1.0_aip/sqrt(fourpi), 0.0_aip, aip) * wingL
        A     = -cmplx(1.0_aip/fourpi, 0.0_aip, aip) * A

        OMP_OFFLOAD target enter data map(to: Binv, wingL, A)
        OMP_OFFLOAD target enter data map(alloc: ag)

        ! Get device pointers
        A_dptr     = get_device_pointer(A    , world%get_device())
        ag_dptr    = get_device_pointer(ag   , world%get_device())
        Binv_dptr  = get_device_pointer(Binv , world%get_device())
        wingL_dptr = get_device_pointer(wingL, world%get_device())

        call _gemm_gpu('n', 'n', nb, 3, nb, -zone, Binv_dptr, &
                        nb, wingL_dptr, nb, zzero, ag_dptr, nb, world)
        call _gemm_gpu('c', 'n', 3, 3, nb, -zone, wingL_dptr, &
                        nb, ag_dptr, nb, zone, A_dptr, 3, world)
        call world%synchronize()

        OMP_OFFLOAD target update from(A, ag)
        OMP_OFFLOAD target exit data map(delete: Binv, wingL, A, ag)

    end subroutine compute_auxiliary_and_A_2d_hermitian


    !> It computes q^T \cdot A \cdot q 
    !> @param[in] A      - A tensor (if GPU it must be loaded)
    !> @param[in] q      - the object q (if GPU it must be loaded)
    !> @param[in] world  - the gpu_world_t handler 
    !> @param[out] qAq   - (if GPU it will also provide the GPU object)
    subroutine compute_qAq(A, q, world, qAq)

        complex(aip), target, allocatable, intent(inout)   :: A(:,:)
        complex(aip), target, allocatable, intent(in)      :: q(:,:)
        type(device_world_t), intent(inout)                :: world
        complex(aip), target, allocatable, intent(out)     :: qAq(:)

        complex(aip), pointer, contiguous :: Aq(:,:)
        integer(i32) :: nr 
        integer(i32) :: i

        type(C_ptr) :: q_dptr, Aq_dptr, A_dptr 

        nr = size(q,2)
        call allocate_device_memory(Aq_dptr, _sizecomplex * product(shape(q)), world%get_device())
        call c_f_pointer(Aq_dptr, Aq, shape(q))
        
        OMP_OFFLOAD target enter data map(to: A)

        ! Get device pointers
        A_dptr   = get_device_pointer(A, world%get_device())
        q_dptr   = get_device_pointer(q, world%get_device())

        call _gemm_gpu('n', 'n', 3, nr, 3, zone, A_dptr, &
                         3, q_dptr, 3, zzero, Aq_dptr, 3, world)
        call world%synchronize()

        allocate(qAq(size(q,2)))

        OMP_OFFLOAD target map(tofrom: qAq) has_device_addr(Aq)
        !$omp teams distribute parallel do private(i)
        do i = 1, size(q,2)
                qAq(i) = dot_product(q(:,i),Aq(:,i))
        end do
        !$omp end teams distribute parallel do
        OMP_OFFLOAD end target

        nullify(Aq)
        call deallocate_device_memory(Aq_dptr, world%get_device())

        OMP_OFFLOAD target exit data map(delete: A)
        deallocate(A)

    end subroutine compute_qAq

    !> It computes ag term for the correction of the body
    !> @param[in] q      - the linear algebra object containing q
    !> @param[in] ag     - ag vector
    !> @param[in] world  - the gpu_world_t handler 
    !> @param[out] qag   - (in GPU if compiled for)
    subroutine compute_inverse_ag(q, ag, world, qag)
        complex(aip), target, allocatable,  intent(inout)   :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)   :: ag(:,:)
        type(device_world_t),       intent(inout)           :: world
        complex(aip), target, allocatable,  intent(inout)   :: qag(:,:)

        integer(i32) :: nr, nb 

        type(C_ptr) :: qag_dptr, ag_dptr, q_dptr
        nr = size(q,2)
        nb = size(ag,1)
        allocate(qag(nr,nb))

        OMP_OFFLOAD target enter data map(to: ag)
        OMP_OFFLOAD target enter data map(alloc: qag)

        ! Get device pointers
        q_dptr    = get_device_pointer(q, world%get_device())
        qag_dptr  = get_device_pointer(qag, world%get_device())
        ag_dptr   = get_device_pointer(ag, world%get_device())

        call _gemm_gpu('t', 't', nr, nb, 3, zone, q_dptr, &
                        3, ag_dptr, nb, zzero, qag_dptr, nr, world)
        call world%synchronize()

        OMP_OFFLOAD target update from(qag)
        OMP_OFFLOAD target exit data map(delete: ag, qag)

        deallocate(ag)

    end subroutine compute_inverse_ag

    !> It computes bg term for the correction of the body
    !> @param[in] q -  qs
    !> @param[in] bg - bg vector
    !> @param[in] world - the gpu_world_t handler 
    !> @param[out] qbg   - (in GPU if compiled for)
    subroutine compute_inverse_bg(q, bg, world, qbg)
        complex(aip), target, allocatable,  intent(inout)  :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)  :: bg(:,:)
        type(device_world_t),               intent(inout)  :: world
        complex(aip), target, allocatable,  intent(out)    :: qbg(:,:)

        integer(i32) :: nr, nb 

        type(C_ptr) :: qbg_dptr, bg_dptr, q_dptr

        nr = size(q,2)
        nb = size(bg,2)
        allocate(qbg(nr,nb))

        OMP_OFFLOAD target enter data map(to: bg)
        OMP_OFFLOAD target enter data map(alloc: qbg)


        ! Get device pointers
        q_dptr    = get_device_pointer(q, world%get_device())
        qbg_dptr  = get_device_pointer(qbg, world%get_device())
        bg_dptr   = get_device_pointer(bg, world%get_device())

        call _gemm_gpu('t', 'n', nr, nb, 3, zone, q_dptr, &
                        3, bg_dptr, 3, zzero, qbg_dptr, nr, world)
        
        call world%synchronize()

        OMP_OFFLOAD target update from(qbg)
        OMP_OFFLOAD target exit data map(delete: bg, qbg)

        deallocate(bg)

    end subroutine compute_inverse_bg

#undef _getrf
#undef _getri
#undef _gemm
#undef _getrf_gpu
#undef _getri_gpu
#undef _gemm_gpu
#undef _sizecomplex

end module idiel_linalg
