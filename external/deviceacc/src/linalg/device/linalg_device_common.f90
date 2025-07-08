! Copyright (C) 2024 DEVICEACC developers
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
!> This file contains unified calls to device accelerated routines for NVIDIA, AMD, and INTEL cards
!> AMD and NVIDIA are supported through MAGMA library, while INTEL cards are supported using INTEL MKL.

!> Module containing unified calls to linear algebra device accelerated routines
module device_linalg_common_interface

    ! The iso_c_binding is included inside the pragmas because its positioning differs in both cases
#if defined(NVIDIAGPU) || defined(AMDGPU)
    use iso_c_binding, only: c_ptr, c_int, c_size_t
    use magma2, only: MagmaLeft, MagmaRight, MagmaLower, MagmaUpper, MagmaNoTrans, &
                      MagmaTrans, MagmaConjTrans, magma_memset_async, &
                      magma_get_zgetri_nb, magma_zgemm, magma_get_cgetri_nb, &
                      magma_zdotu, magma_cdotu, magma_zdotc, magma_cdotc, &
                      magma_cgerc, magma_caxpy, magmablas_cgeadd2, magma_cgetri_gpu, &
                      magmablas_zgeadd2, magma_zgemm_batched_strided, magma_zhemv, magma_chemv, &
                      magma_chemm, magma_ccopyvector_async, magma_zgeru, magma_zcopyvector_async, &
                      magma_zhemm, magma_cgetrf_gpu, magma_chemv, magma_zgerc, &
                      magma_zgetrf_gpu, magma_cgeru, magma_cgemm_batched_strided, magma_zgetri_gpu, &
                      magma_zaxpy, magma_cgemm
#endif
#if defined(INTELGPU)
    use iso_c_binding, only: c_ptr, c_int, c_size_t, c_f_pointer
    use onemkl_lapack_omp_offload_lp64, only: cgetrf, zgetrf, cgetri, zgetri
    use onemkl_blas_omp_offload_lp64,   only: cgemm, zgemm, chemv, zhemv, chemm, zhemm, &
                                              cgemm_batch_strided, zgemm_batch_strided, &
                                              cdotu, zdotu, cdotc, zdotc, caxpy, zaxpy, &
                                              cgerc, zgerc, cgeru, zgeru, mkl_comatadd_batch_strided, &
                                              mkl_zomatadd_batch_strided, ccopy, zcopy
#endif
    use iso_fortran_env,  only: i32=>int32, r32=>real32, r64=>real64
    use m_memory_device,  only: bytes_single_complex, bytes_double_complex
    use m_device_world_t, only: device_world_t

    implicit none

    private 
    public :: cgemm_gpu,  cdotc_gpu, cdotu_gpu, cgetrf_gpu, cgetri_gpu, get_cgetri_nb_gpu, &
              zgemm_gpu,  zdotc_gpu, zdotu_gpu, zgetrf_gpu, zgetri_gpu, get_zgetri_nb_gpu, &
              caxpy_gpu,  zaxpy_gpu, cgemm_batched_gpu, zgemm_batched_gpu, &
              cgerc_gpu,  cgeru_gpu, zgerc_gpu, zgeru_gpu, &
              chemv_gpu,  zhemv_gpu, chemm_gpu, zhemm_gpu, &
              cgeadd_gpu, zgeadd_gpu, &
              ccopy_gpu,  zcopy_gpu, &
              csetzero_gpu, zsetzero_gpu

contains

    !!!!!!!!!!!!!!!  SINGLE PRECISION !!!!!!!!!!!!!!

    !> Complex single precision LU decomposition.
    !> @param[in] m - The number of rows of the matrix A
    !> @param[in] n - The number of columns of the matrix A
    !> @param[in,out] dA - C-pointer to A matrix (entry) in the device. On exit, is a C pointer to the factors L and U from the factorization.
    !> @param[in] lda - The leading dimension of the array A.
    !> @param[in,out] ipiv - The pivot indices; for 1 <= i <= min(m,n), row i of the matrix was interchanged with row ipiv(i). Host array.
    !> @param[out] info - execution information. 0 success; < 0:  if INFO = -i, the i-th argument had an illegal value;  > 0:  if INFO = i, U(i,i) is exactly zero.
    !> @param[in,out] world - the device-host handler
    subroutine cgetrf_gpu(m, n, dA, lda, ipiv, info, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_cgetrf_gpu(m, n, dA, lda, ipiv, info)
#endif
#if defined(INTELGPU)
        complex(r32), pointer :: A(:,:)

        !We associate A with the device pointer dA, pointing to GPU memory.
        !This association allows the Fortran code to access and manipulate the data
        !on the GPU without explicitly transferring it back and forth between the CPU and GPU.
        call c_f_pointer(dA, A, [lda,n])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch is_device_ptr(A)
        call cgetrf(m, n, A, lda, ipiv, info)
        !$omp end target data

        nullify(A)
#endif
    end subroutine cgetrf_gpu

    !> Computes the inverse of a single precision complex matrix using its LU factorization.
    !> @param[in] n - The order of the matrix A
    !> @param[in,out] dA - On entry, the C-pointer to factors L and U from the factorization of matrix A. On exit, if successful, contains the address to the inverse of A. Host C-pointer to device memory..
    !> @param[in] ipiv - pivot indices. Host array.
    !> @param[in,out] dwork - Host C-pointer to device memory to work space.
    !> @param[in,out] lwork - size of the workspace
    !> @param[in,out] world - the device-host handler
    subroutine cgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(c_ptr),  value                             :: dwork
        integer(i32), intent(in)                        :: lwork
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_cgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info)
#endif
#if defined(INTELGPU)

        complex(r32), pointer :: A(:,:), work(:)

        call c_f_pointer(dA, A, [lda,n])
        call c_f_pointer(dwork, work, [lwork])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch is_device_ptr(A)
        call cgetri(n, A, lda, ipiv, work, lwork, info)
        !$omp end target data

        nullify(A, work)
#endif

    end subroutine cgetri_gpu

    !> Provides the appropiate block size for the workspace of the single precision complex matrix inversion using LU decomposition
    !> The optimal worksize is then obtained by multiplying n by the result of this function
    !> @param[in] n - the order of the A matrix
    function get_cgetri_nb_gpu(n) result(nblock)
        integer(i32), intent(in) :: n

        integer :: nblock

#if defined(NVIDIAGPU) || defined(AMDGPU)
        nblock = magma_get_cgetri_nb(n)
#endif
#if defined (INTELGPU)
        nblock = 64_i32
#endif

    end function get_cgetri_nb_gpu

    !> Complex single precision matrix-matrix product (C = alpha * op(A) * op(B) + beta * C)
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      -  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] dA     - C-pointer to the A matrix. Host C-pointer to device memory
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] dB     - C-pointer to the B matrix. Host C-pointer to device memory
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] dC - C-pointer to the C matrix. Host C-pointer to device memory
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in,out] world - device-host handler.
    subroutine cgemm_gpu(transa, transb, m, n, k, alpha, dA, lda, dB, ldb, beta, dC, ldc, world)
        character, intent(in)         :: transa
        character, intent(in)         :: transb
        integer(i32), intent(in)      :: m
        integer(i32), intent(in)      :: n
        integer(i32), intent(in)      :: k
        complex(r32), intent(in)      :: alpha
        type(c_ptr),  value           :: dA
        integer(i32), intent(in)      :: lda
        type(c_ptr),  value           :: dB
        integer(i32), intent(in)      :: ldb
        complex(r32), intent(in)      :: beta
        type(c_ptr),  value           :: dC
        integer(i32), intent(in)      :: ldc
        type(device_world_t), intent(in) :: world

        integer(c_int) :: opA, opB

#if defined(NVIDIAGPU) || defined(AMDGPU)

        if (transa == 'N' .or. transa == 'n') opA = MagmaNoTrans
        if (transa == 'T' .or. transa == 't') opA = MagmaTrans
        if (transa == 'C' .or. transa == 'c') opA = MagmaConjTrans

        if (transb == 'N' .or. transb == 'n') opB = MagmaNoTrans
        if (transb == 'T' .or. transb == 't') opB = MagmaTrans
        if (transb == 'C' .or. transb == 'c') opB = MagmaConjTrans

        call magma_cgemm(opA, opB, m, n, k, alpha, dA, lda, dB, ldb, beta, dC, ldc, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)
        integer(i32) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(dA, A, [lda,ka])
        call c_f_pointer(dB, B, [ldb,kb])
        call c_f_pointer(dC, C, [ldc,n])

        !$omp dispatch is_device_ptr(A,B,C)
        call cgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

        nullify(A, B, C)

#endif
    end subroutine cgemm_gpu

    !> Complex single precision matrix-matrix batched product (C[ibatch] = alpha * op(A[ibatch]) * op(B[ibatch]) + beta * C[ibatch])
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      -  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] dA     - C-pointer to the A matrices. Host C-pointer to device memory
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] dB     - C-pointer to the B matrices. Host C-pointer to device memory
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] dC - C-pointer to the C matrices. Host C-pointer to device memory
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in] batchcount - number of matrices in the batch
    !> @param[in,out] world - device-host handler.
    !> @param[in,optional] stridea - If present it indicates a different stride for A matrix than the matrix size.
    !> @param[in,optional] strideb - If present it indicates a different stride for B matrix than the matrix size.
    !> @param[in,optional] stridec - If present it indicates a different stride for C matrix than the matrix size.
    subroutine cgemm_batched_gpu(transa, transb, m, n, k, alpha, dA, lda, dB, ldb, beta, dC, ldc, batchcount, world, stridea, &
                                 strideb, stridec)
        character, intent(in)              :: transa
        character, intent(in)              :: transb
        integer(i32), intent(in)           :: m
        integer(i32), intent(in)           :: n
        integer(i32), intent(in)           :: k
        complex(r32), intent(in)           :: alpha
        type(c_ptr),  value                :: dA
        integer(i32), intent(in)           :: lda
        type(c_ptr),  value                :: dB
        integer(i32), intent(in)           :: ldb
        complex(r32), intent(in)           :: beta
        type(c_ptr),  value                :: dC
        integer(i32), intent(in)           :: ldc
        integer(i32), intent(in)           :: batchcount
        type(device_world_t), intent(in)   :: world
        integer(i32), intent(in), optional :: stridea
        integer(i32), intent(in), optional :: strideb
        integer(i32), intent(in), optional :: stridec

        integer(c_int) :: opA, opB
        integer(c_int) :: ka, kb

        integer(i32) :: stridea_local, strideb_local, stridec_local
        type(c_ptr), target :: dA_array(batchcount), dB_array(batchcount), dC_array(batchcount)
        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        stridea_local = lda*ka
        strideb_local = ldb*kb
        stridec_local = ldc*n
        if (present(stridea)) stridea_local = stridea
        if (present(strideb)) strideb_local = strideb
        if (present(stridec)) stridec_local = stridec

#if defined(NVIDIAGPU) || defined(AMDGPU)
        
        if (transa == 'N' .or. transa == 'n') opA = MagmaNoTrans
        if (transa == 'T' .or. transa == 't') opA = MagmaTrans
        if (transa == 'C' .or. transa == 'c') opA = MagmaConjTrans

        if (transb == 'N' .or. transb == 'n') opB = MagmaNoTrans
        if (transb == 'T' .or. transb == 't') opB = MagmaTrans
        if (transb == 'C' .or. transb == 'c') opB = MagmaConjTrans

        call magma_cgemm_batched_strided(opA, opB, m, n, k, alpha, dA, lda, stridea_local, &
                                         dB, ldb, strideb_local, &
                                         beta, dC, ldc, stridec_local, &
                                         batchcount, world%get_queue())

#endif
#if defined(INTELGPU)

        call c_f_pointer(dA, A, [lda,ka])
        call c_f_pointer(dB, B, [ldb,kb])
        call c_f_pointer(dC, C, [ldc,n])

        !$omp dispatch is_device_ptr(A,B,C)
        call cgemm_batch_strided(transa, transb, m, n, k, alpha, A, lda, stridea_local, B, ldb, &
                                   strideb_local, beta, C, ldc, stridec_local, batchcount)

        nullify(A, B, C)

#endif
    end subroutine cgemm_batched_gpu

    !> Complex single precision dot product (unconjugated) of vectors x and y; \( x^T y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Host C-pointer to device memory
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Host C-pointer to device memory
    !> @param[in] incy - Stride between consecutive elements of dy
    !> @param[in,out] world - device-host handler.
    complex(r32) function cdotu_gpu(n, dx, incx, dy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: dx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: dy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world
 
#if defined(NVIDIAGPU) || defined(AMDGPU)
        cdotu_gpu = magma_cdotu(n, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)
        complex(r32), pointer :: x(:), y(:)
        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        !$omp dispatch is_device_ptr(x,y)
        cdotu_gpu = cdotu(n, x, incx, y, incy)
        nullify(x, y)
#endif 
    end function cdotu_gpu

    !> Complex single precision dot product (conjugated) of vectors x and y; \( x^H y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Host C-pointer to device memory
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Host C-pointer to device memory
    !> @param[in] incy - Stride between consecutive elements of dy
    !> @param[in,out] world - device-host handler.
    complex(r32) function cdotc_gpu(n, dx, incx, dy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: dx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: dy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world
        
#if defined(NVIDIAGPU) || defined(AMDGPU)
        cdotc_gpu = magma_cdotc(n, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)
        complex(r32), pointer :: x(:), y(:)
        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        
        !$omp dispatch is_device_ptr(x,y)
        cdotc_gpu = cdotc(n, x, incx, y, incy)
        
        nullify(x, y)
#endif 
    end function cdotc_gpu

    !> Complex single precision constant times a vector plus a vector; \( y = \alpha x + y \). 
    !> @param[in]	n	- Number of elements in vectors x and y. n >= 0.
    !> @param[in]	alpha	- Scalar \( \alpha \)
    !> @param[in]	dx	- Host C-pointer to device memory to x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]	incx	- Stride between consecutive elements of dx. incx != 0.
    !> @param[in,out]	dy	- Host C-pointer to device memory to y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]	incy	- Stride between consecutive elements of dy. incy != 0.
    !> @param[in,out]   world	- the device-host handler
    subroutine caxpy_gpu(n, alpha, dx, incx, dy, incy, world)
        integer(i32), intent(in)                        :: n
        complex(r32), intent(in)                        :: alpha
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_caxpy(n, alpha, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r32), pointer :: x(:), y(:)

        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])

        !$omp dispatch is_device_ptr(x,y)
        call caxpy(n, alpha, x, incx, y, incy)

        nullify(x,y)

#endif 

    end subroutine caxpy_gpu
            
    !> Perform rank-1 update, \( A = \alpha x y^H + A \) for complex double precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] dx -    Host C-pointer to device memory to vector x. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] dy - Host C-pointer to device memory to vector y. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y. Must not be zero.
    !> @param[in,out] dA - Host C-pointer to device memory to A. On entry, the m-by-n matrix A. On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldda >= max(1, m))
    !> @param[in,out]   world - the device-host handler 
    subroutine cgerc_gpu(m, n, alpha, dx, incx, dy, incy, dA, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r32), intent(in)                        :: alpha
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_cgerc(m, n, alpha, dx, incx, dy, incy, dA, lda, world%get_queue())
#endif

#if defined(INTELGPU)
        
        complex(r32), pointer :: x(:), y(:), a(:,:)

        call c_f_pointer(dx, x, [1 + (m-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        call c_f_pointer(dA, A, [lda,n])
        !$omp dispatch is_device_ptr(x,y,A)
        call cgerc(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)
#endif

    end subroutine cgerc_gpu

    !> Perform rank-1 update, \( A = \alpha x y^T + A \) for complex single precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] dx -    Host C-pointer to device memory to vector x. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] dy - Host C-pointer to device memory to vector y. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y. Must not be zero.
    !> @param[in,out] dA - Host C-pointer to device memory to A. On entry, the m-by-n matrix A. On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldda >= max(1, m))
    !> @param[in,out]   world - the device-host handler 
    subroutine cgeru_gpu(m, n, alpha, dx, incx, dy, incy, dA, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r32), intent(in)                        :: alpha
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_cgeru(m, n, alpha, dx, incx, dy, incy, dA, lda, world%get_queue())
#endif
        
#if defined(INTELGPU)
        
        complex(r32), pointer :: x(:), y(:), A(:,:)

        call c_f_pointer(dx, x, [1 + (m-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        call c_f_pointer(dA, A, [lda,n])
        !$omp dispatch is_device_ptr(x,y,A)
        call cgeru(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)
#endif

    end subroutine cgeru_gpu

    !> Perform Hermitian matrix-vector product in single precision.
    !> \( y = \alpha A x + \beta y \) (side == MagmaLeft), or
    !> \( y = \alpha A v + \beta y \) (side == MagmaRight),
    !> where \( A \) is Hermitian.
    !>
    !> @param[in] uplo  - specifies  whether  the  upper ('U' or 'u')  or  lower ('L' or 'l') triangular  part  of  the  hermitian matrix is used
    !> @param[in] n     - the order of the A matrix
    !> @param[in] alpha - the \(\alpha\) scaling factor
    !> @param[in] dA    - Host C-pointer to device memory to matrix A
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] dx    - Host C-pointer to device memory to vector x
    !> @param[in] incx  - Specifies the increment for the elements of x.
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] dy    - Host C-pointer to device memory to vector y
    !> @param[in] incy  - Specifies the increment for the elements of y.
    !> @param[in] world - the device-host handler
    subroutine chemv_gpu(uplo, n, alpha, dA, lda, dx, incx, beta, dy, incy, world)
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: n
        complex(r32), intent(in)            :: alpha
        type(c_ptr),  value                 :: dA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: dx
        integer(i32), intent(in)            :: incx
        complex(r32), intent(in)            :: beta
        type(c_ptr),  value                 :: dy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world

        integer(i32) :: magma_side, magma_uplo
        complex(r32), pointer :: A(:,:), x(:), y(:)

#if defined(NVIDIAGPU) || defined(AMDGPU)

        magma_uplo = merge(MagmaLower, MagmaUpper, uplo == 'L' .or. uplo == 'l')

        call magma_chemv(magma_uplo, n, alpha, dA, lda, dx, incx, beta, dy, incy, world%get_queue())

#endif

#if defined(INTELGPU)

        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        call c_f_pointer(dA, A, [lda,n])

        !$omp dispatch is_device_ptr(A,x,y)
        call chemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)

        nullify(A,x,y)
#endif

    end subroutine chemv_gpu

    !> Perform Hermitian matrix-matrix product in single precision.
    !> \( C = \alpha A B + \beta C \) (side == MagmaLeft), or
    !> \( C = \alpha B A + \beta C \) (side == MagmaRight),
    !> where \( A \) is Hermitian.
    !>
    !> @param[in] side  - Whether A is on the left ('L' or 'l') or right ('R', 'r').
    !> @param[in] uplo  - specifies  whether  the  upper ('U' or 'u')  or  lower ('L' or 'l') triangular  part  of  the  hermitian matrix is used
    !> @param[in] m     - the number of rows of C
    !> @param[in] n     - the number of columns of C
    !> @param[in] alpha - the \(\alpha\) scaling factor
    !> @param[in] dA    - Host C-pointer to device memory to matrix A
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] dB    - Host C-pointer to device memory to matrix B
    !> @param[in] ldb   - leading dimension of the matrix B
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] dC    - Host C-pointer to device memory to matrix C
    !> @param[in] ldc   - leading dimension of the matrix C
    !> @param[in] world - the device-host handler
    subroutine chemm_gpu(side, uplo, m, n, alpha, dA, lda, dB, ldb, beta, dC, ldc, world)
        character, intent(in)               :: side
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r32), intent(in)            :: alpha
        type(c_ptr),  value                 :: dA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: dB
        integer(i32), intent(in)            :: ldb
        complex(r32), intent(in)            :: beta
        type(c_ptr),  value                 :: dC
        integer(i32), intent(in)            :: ldc
        type(device_world_t), intent(inout) :: world

        integer(i32) :: ka
        integer(i32) :: magma_side, magma_uplo
        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)

#if defined(NVIDIAGPU) || defined(AMDGPU)

        magma_side = merge(MagmaLeft, MagmaRight, side == 'L' .or. side == 'l')
        magma_uplo = merge(MagmaLower, MagmaUpper, uplo == 'L' .or. uplo == 'l')

        call magma_chemm(magma_side, magma_uplo, m, n, alpha, dA, lda, dB, ldb, beta, dC, ldc, world%get_queue())

#endif

#if defined(INTELGPU)

        ka = merge(m,n,side == 'L' .or. side == 'l')

        call c_f_pointer(dA, A, [lda,ka])
        call c_f_pointer(dB, B, [ldb,n])
        call c_f_pointer(dC, C, [ldc,n])

        !$omp dispatch is_device_ptr(A,B,C)
        call chemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)

        nullify(A,B,C)
#endif

    end subroutine chemm_gpu 

    !> Computes dB = alpha*dA + beta*dB. Single precision
    !> @param[in] m     - The number of rows of the matrix A
    !> @param[in] n     - The number of columns of the matrix A
    !> @param[in] alpha - alpha scaling factor.
    !> @param[in] dA    - Host C-pointer to A matrix.
    !> @param[in] lda   - The leading dimension of the array A.
    !> @param[in] beta  - beta scaling factor.
    !> @param[in,out] dB    - Host C-pointer to B matrix. On exit it have the result.
    !> @param[in] ldb   - The leading dimension of the array B.
    !> @param[in,out] world - the device-host handler. 
    subroutine cgeadd_gpu(m, n, alpha, dA, lda, beta, dB, ldb, world)
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r32), intent(in)            :: alpha
        type(c_ptr),  value                 :: dA
        integer(i32), intent(in)            :: lda
        complex(r32), intent(in)            :: beta
        type(c_ptr),  value                 :: dB
        integer(i32), intent(in)            :: ldb
        type(device_world_t), intent(inout) :: world

#if defined(INTELGPU)
        complex(r32), pointer :: A(:,:), B(:,:)
        call c_f_pointer(dA, A, int([lda,n],c_size_t))
        call c_f_pointer(dB, B, int([ldb,n],c_size_t))
        !$omp dispatch is_device_ptr(A,B)
        call mkl_comatadd_batch_strided('c','n','n',m,n,alpha,A,lda,lda*m,beta,B,ldb,ldb*m,B,ldb,ldb*m,1)
        nullify(A,B)
#endif

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magmablas_cgeadd2(m, n, alpha, dA, lda, beta, dB, ldb, world%get_queue())
#endif

    end subroutine cgeadd_gpu

    !> Copies x contents to y.
    !> @param[in]       n           - Number of elements in vectors x and y. n >= 0.
    !> @param[in]       dx          - Host pointer to device x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]       incx        - Stride between consecutive elements of dx.
    !> @param[in,out]   dy          - Host pointer to device y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]       incy        - Stride between consecutive elements of dy.
    !> @param[in,out]   world       - the device-host handler
    subroutine ccopy_gpu(n, dx, incx, dy, incy, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_ccopyvector_async(n, dx, incx, dy, incy, world%get_queue())
#endif

#if defined(INTELGPU)
        complex(r32), pointer :: x(:), y(:)
        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        !$omp dispatch is_device_ptr(x,y)
        call ccopy(n, x, incx, y, incy)
        nullify(x,y)
#endif

    end subroutine ccopy_gpu

    !> Sets an array to zero
    !> @param[in,out] dA    - Host pointer to device A array.
    !> @param[in] size  - the size of the array
    !> @param[in,out] world - the device-host handler. 
    subroutine csetzero_gpu(dA, size, world)
        type(c_ptr),  value                 :: dA
        integer(c_size_t), intent(in)       :: size
        type(device_world_t), intent(inout) :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        integer(i32) :: ierr
        ierr = magma_memset_async(dA, 0_c_int, size*bytes_single_complex, world%get_queue())
#endif

#if defined(INTELGPU)
        complex(r32), pointer   :: a(:)
        complex(r32), parameter :: zzero = cmplx(0.0_r64,0.0_r64,kind=r64)
        integer(c_size_t) :: i
        call c_f_pointer(dA, a, [size])
        !$omp target has_device_addr(a)
        !$omp teams distribute parallel do simd
        do i = 1, size
            a(i) = zzero
        end do
        !$omp end teams distribute parallel do simd
        !$omp end target
        nullify(a)
#endif

    end subroutine csetzero_gpu

    !!!!!!!!!!!!!!!  DOUBLE PRECISION !!!!!!!!!!!!!!
    
    !> Complex double precision LU decomposition.
    !> @param[in] m - The number of rows of the matrix A
    !> @param[in] n - The number of columns of the matrix A
    !> @param[in,out] dA - C pointer to A matrix (entry) in the device. On exit, is a C pointer to the factors L and U from the factorization.
    !> @param[in] lda - The leading dimension of the array A.
    !> @param[in,out] ipiv - The pivot indices; for 1 <= i <= min(m,n), row i of the matrix was interchanged with row ipiv(i). Host array.
    !> @param[out] info - execution information. 0 success; < 0:  if INFO = -i, the i-th argument had an illegal value;  > 0:  if INFO = i, U(i,i) is exactly zero.
    !> @param[in,out] world - the device-host handler
    subroutine zgetrf_gpu(m, n, dA, lda, ipiv, info, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zgetrf_gpu(m, n, dA, lda, ipiv, info)
#endif
#if defined(INTELGPU)
        complex(r64), pointer :: A(:,:)
        ! We associate A with the device pointer dA, pointing to GPU memory.
        ! This association allows the Fortran code to access and manipulate the data
        ! on the GPU without explicitly transferring it back and forth between the CPU and GPU.
        call c_f_pointer(dA, A, [lda,n])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch is_device_ptr(A)
        call zgetrf(m, n, A, lda, ipiv, info)
        !$omp end target data

        nullify(A)
#endif
    end subroutine zgetrf_gpu

    !> Computes the inverse of a double precision complex matrix using its LU factorization.
    !> @param[in] n - The order of the matrix A
    !> @param[in,out] dA - On entry, the C-pointer to factors L and U from the factorization of matrix A. On exit, if successful, contains the address to the inverse of A. Host C-pointer to device memory..
    !> @param[in] ipiv - pivot indices. Host array.
    !> @param[in,out] dwork - Host C-pointer to device memory to work space.
    !> @param[in,out] lwork - size of the workspace
    !> @param[in,out] world - the device-host handler
    subroutine zgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(c_ptr),  value                             :: dwork
        integer(i32), intent(in)                        :: lwork
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info)
#endif
#if defined(INTELGPU)
        complex(r64), pointer :: A(:,:), work(:)

        call c_f_pointer(dA, A, [lda,n])
        call c_f_pointer(dwork, work, [lwork])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch is_device_ptr(A)
        call zgetri(n, A, lda, ipiv, work, lwork, info)
        !$omp end target data

        nullify(A, work)
#endif

    end subroutine zgetri_gpu

    !> Provides the appropiate block size for the workspace of the double precision complex matrix inversion using LU decomposition
    !> The optimal worksize is then obtained by multiplying n by the result of this function
    !> @param[in] n - the order of the A matrix
    function get_zgetri_nb_gpu(n) result(nblock)
        integer(i32), intent(in) :: n

        integer :: nblock

#if defined(NVIDIAGPU) || defined(AMDGPU)
        nblock = magma_get_zgetri_nb(n)
#endif
#if defined (INTELGPU)
        nblock = 64_i32
#endif

    end function get_zgetri_nb_gpu

    !> Complex double precision matrix-matrix product (C = alpha * op(A) * op(B) + beta * C)
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      -  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] da     - C-pointer to the A matrix. Host C-pointer to device memory
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] db     - C-pointer to the B matrix. Host C-pointer to device memory
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] dc - C-pointer to the C matrix. Host C-pointer to device memory
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in,out] world - device-host handler.
    subroutine zgemm_gpu(transa, transb, m, n, k, alpha, dA, lda, dB, ldb, beta, dC, ldc, world)
        character, intent(in)            :: transa
        character, intent(in)            :: transb
        integer(i32), intent(in)         :: m
        integer(i32), intent(in)         :: n
        integer(i32), intent(in)         :: k
        complex(r64), intent(in)         :: alpha
        type(c_ptr),  value              :: dA
        integer(i32), intent(in)         :: lda
        type(c_ptr),  value              :: dB
        integer(i32), intent(in)         :: ldb
        complex(r64), intent(in)         :: beta
        type(c_ptr),  value              :: dC
        integer(i32), intent(in)         :: ldc
        type(device_world_t), intent(in) :: world

        integer(c_int) :: opA, opB

#if defined(NVIDIAGPU) || defined(AMDGPU)

        if (transa == 'N' .or. transa == 'n') opA = MagmaNoTrans
        if (transa == 'T' .or. transa == 't') opA = MagmaTrans
        if (transa == 'C' .or. transa == 'c') opA = MagmaConjTrans

        if (transb == 'N' .or. transb == 'n') opB = MagmaNoTrans
        if (transb == 'T' .or. transb == 't') opB = MagmaTrans
        if (transb == 'C' .or. transb == 'c') opB = MagmaConjTrans

        call magma_zgemm(opA, opB, m, n, k, alpha, dA, lda, dB, ldb, beta, dC, ldc, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)
        integer(i32) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(dA, A, [lda,ka])
        call c_f_pointer(dB, B, [ldb,kb])
        call c_f_pointer(dC, C, [ldc,n])

        !$omp dispatch is_device_ptr(A,B,C) 
        call zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

        nullify(A, B, C)

#endif
    end subroutine zgemm_gpu

    !> Complex double precision matrix-matrix batched product (C[ibatch] = alpha * op(A[ibatch]) * op(B[ibatch]) + beta * C[ibatch])
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      -  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] dA     - C-pointer to the A matrices. Host C-pointer to device memory
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] dB     - C-pointer to the B matrices. Host C-pointer to device memory
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] dC - C-pointer to the C matrices. Host C-pointer to device memory
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in] batchcount - number of matrices in the batch
    !> @param[in,out] world - device-host handler.
    !> @param[in,optional] stridea - If present it indicates a different stride for A matrix than the matrix size.
    !> @param[in,optional] strideb - If present it indicates a different stride for B matrix than the matrix size.
    !> @param[in,optional] stridec - If present it indicates a different stride for C matrix than the matrix size.
    subroutine zgemm_batched_gpu(transa, transb, m, n, k, alpha, dA, lda, dB, ldb, beta, dC, ldc, batchcount, world, stridea, &
                                 strideb, stridec)
        character, intent(in)              :: transa
        character, intent(in)              :: transb
        integer(i32), intent(in)           :: m
        integer(i32), intent(in)           :: n
        integer(i32), intent(in)           :: k
        complex(r64), intent(in)           :: alpha
        type(c_ptr),  value                :: dA
        integer(i32), intent(in)           :: lda
        type(c_ptr),  value                :: dB
        integer(i32), intent(in)           :: ldb
        complex(r64), intent(in)           :: beta
        type(c_ptr),  value                :: dC
        integer(i32), intent(in)           :: ldc
        integer(i32), intent(in)           :: batchcount
        type(device_world_t), intent(in)   :: world
        integer(i32), intent(in), optional :: stridea
        integer(i32), intent(in), optional :: strideb
        integer(i32), intent(in), optional :: stridec

        integer(c_int) :: opA, opB
        integer(c_int) :: ka, kb

        integer(i32) :: stridea_local, strideb_local, stridec_local
        type(c_ptr), target :: dA_array(batchcount), dB_array(batchcount), dC_array(batchcount)
        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        stridea_local = lda*ka
        strideb_local = ldb*kb
        stridec_local = ldc*n
        if (present(stridea)) stridea_local = stridea
        if (present(strideb)) strideb_local = strideb
        if (present(stridec)) stridec_local = stridec

#if defined(NVIDIAGPU) || defined(AMDGPU)
        
        if (transa == 'N' .or. transa == 'n') opA = MagmaNoTrans
        if (transa == 'T' .or. transa == 't') opA = MagmaTrans
        if (transa == 'C' .or. transa == 'c') opA = MagmaConjTrans

        if (transb == 'N' .or. transb == 'n') opB = MagmaNoTrans
        if (transb == 'T' .or. transb == 't') opB = MagmaTrans
        if (transb == 'C' .or. transb == 'c') opB = MagmaConjTrans

        call magma_zgemm_batched_strided(opA, opB, m, n, k, alpha, dA, lda, stridea_local, &
                                         dB, ldb, strideb_local, &
                                         beta, dC, ldc, stridec_local, &
                                         batchcount, world%get_queue())


#endif
#if defined(INTELGPU)

        call c_f_pointer(dA, A, [lda,ka])
        call c_f_pointer(dB, B, [ldb,kb])
        call c_f_pointer(dC, C, [ldc,n])

        !$omp dispatch is_device_ptr(A,B,C)
        call zgemm_batch_strided(transa, transb, m, n, k, alpha, A, lda, stridea_local, B, ldb, &
                                   strideb_local, beta, C, ldc, stridec_local, batchcount)

        nullify(A,B,C)

#endif
    end subroutine zgemm_batched_gpu

    !> Complex double precision dot product (unconjugated) of vectors x and y; \( x^T y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Host C-pointer to device memory
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Host C-pointer to device memory
    !> @param[in] incy - Stride between consecutive elements of dy
    !> @param[in,out] world - device-host handler.
    complex(r64) function zdotu_gpu(n, dx, incx, dy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: dx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: dy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world
 
#if defined(NVIDIAGPU) || defined(AMDGPU)
        zdotu_gpu = magma_zdotu(n, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)
        complex(r64), pointer :: x(:), y(:)
        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        !$omp dispatch is_device_ptr(x,y)
        zdotu_gpu = zdotu(n, x, incx, y, incy)
        nullify(x, y)
#endif 
    end function zdotu_gpu

    !> Complex double precision dot product (conjugated) of vectors x and y; \( x^H y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Host C-pointer to device memory
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Host C-pointer to device memory
    !> @param[in] incy - Stride between consecutive elements of dy
    !> @param[in,out] world - device-host handler.
    complex(r64) function zdotc_gpu(n, dx, incx, dy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: dx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: dy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world
        
#if defined(NVIDIAGPU) || defined(AMDGPU)
        zdotc_gpu = magma_zdotc(n, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)
        complex(r64), pointer :: x(:), y(:)
        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        
        !$omp dispatch is_device_ptr(x,y)
        zdotc_gpu = zdotc(n, x, incx, y, incy)
        
        nullify(x, y)
#endif 
    end function zdotc_gpu

    !> Complex double precision constant times a vector plus a vector; \( y = \alpha x + y \). 
    !> @param[in]	n	- Number of elements in vectors x and y. n >= 0.
    !> @param[in]	alpha	- Scalar \( \alpha \)
    !> @param[in]	dx	- Host C-pointer to device memory to x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]	incx	- Stride between consecutive elements of dx. incx != 0.
    !> @param[in,out]	dy	- Host C-pointer to device memory to y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]	incy	- Stride between consecutive elements of dy. incy != 0.
    !> @param[in,out]   world	- the device-host handler
    subroutine zaxpy_gpu(n, alpha, dx, incx, dy, incy, world)
        integer(i32), intent(in)                        :: n
        complex(r64), intent(in)                        :: alpha
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zaxpy(n, alpha, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r64), pointer :: x(:), y(:)

        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])

        !$omp dispatch is_device_ptr(x,y)
        call zaxpy(n, alpha, x, incx, y, incy)

        nullify(x,y)

#endif 

    end subroutine zaxpy_gpu

    !> Perform rank-1 update, \( A = \alpha x y^H + A \) for complex double precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] dx -    Host C-pointer to device memory to vector x. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] dy - Host C-pointer to device memory to vector y. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y. Must not be zero.
    !> @param[in,out] dA - Host C-pointer to device memory to A. On entry, the m-by-n matrix A. On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldda >= max(1, m))
    !> @param[in,out]   world - the device-host handler 
    subroutine zgerc_gpu(m, n, alpha, dx, incx, dy, incy, dA, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r64), intent(in)                        :: alpha
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zgerc(m, n, alpha, dx, incx, dy, incy, dA, lda, world%get_queue())
#endif

#if defined(INTELGPU)
        
        complex(r64), pointer :: x(:), y(:), A(:,:)

        call c_f_pointer(dx, x, [1 + (m-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        call c_f_pointer(dA, A, [lda,n])
        !$omp dispatch is_device_ptr(x,y,A)
        call zgerc(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)
#endif

    end subroutine zgerc_gpu

    !> Perform rank-1 update, \( A = \alpha x y^T + A \) for complex double precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] dx -    Host pointer to device vector x. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] dy - Host pointer to device vector y. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y. Must not be zero.
    !> @param[in,out] dA - Host pointer to device A. On entry, the m-by-n matrix A. On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldda >= max(1, m))
    !> @param[in,out]   world - the device-host handler 
    subroutine zgeru_gpu(m, n, alpha, dx, incx, dy, incy, da, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r64), intent(in)                        :: alpha
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zgeru(m, n, alpha, dx, incx, dy, incy, dA, lda, world%get_queue())
#endif
        
#if defined(INTELGPU)
        
        complex(r64), pointer :: x(:), y(:), A(:,:)

        call c_f_pointer(dx, x, [1 + (m-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        call c_f_pointer(dA, A, [lda,n])
        !$omp dispatch is_device_ptr(x,y,A)
        call zgeru(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)
#endif

    end subroutine zgeru_gpu

    !> Perform Hermitian matrix-vector product in single precision.
    !> \( y = \alpha A x + \beta y \) (side == MagmaLeft), or
    !> \( y = \alpha A v + \beta y \) (side == MagmaRight),
    !> where \( A \) is Hermitian.
    !>
    !> @param[in] uplo  - specifies  whether  the  upper ('U' or 'u')  or  lower ('L' or 'l') triangular  part  of  the  hermitian matrix is used
    !> @param[in] n     - the order of the A matrix
    !> @param[in] alpha - the \(\alpha\) scaling factor
    !> @param[in] da    - Host C-pointer to device memory to matrix A
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] dx    - Host C-pointer to device memory to vector x
    !> @param[in] incx  - Specifies the increment for the elements of x.
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] dy    - Host C-pointer to device memory to vector y
    !> @param[in] incy  - Specifies the increment for the elements of y.
    !> @param[in] world - the device-host handler
    subroutine zhemv_gpu(uplo, n, alpha, da, lda, dx, incx, beta, dy, incy, world)
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: n
        complex(r64), intent(in)            :: alpha
        type(c_ptr),  value                 :: da
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: dx
        integer(i32), intent(in)            :: incx
        complex(r64), intent(in)            :: beta
        type(c_ptr),  value                 :: dy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world

        integer(i32) :: magma_side, magma_uplo
        complex(r64), pointer :: a(:,:), x(:), y(:)

#if defined(NVIDIAGPU) || defined(AMDGPU)

        magma_uplo = merge(MagmaLower, MagmaUpper, uplo == 'L' .or. uplo == 'l')

        call magma_zhemv(magma_uplo, n, alpha, da, lda, dx, incx, beta, dy, incy, world%get_queue())

#endif

#if defined(INTELGPU)

        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        call c_f_pointer(da, a, [lda,n])

        !$omp dispatch is_device_ptr(a,x,y)
        call zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)

        nullify(a,x,y)
#endif

    end subroutine zhemv_gpu

    !> Perform Hermitian matrix-matrix product in single precision.
    !> \( C = \alpha A B + \beta C \) (side == MagmaLeft), or
    !> \( C = \alpha B A + \beta C \) (side == MagmaRight),
    !> where \( A \) is Hermitian.
    !>
    !> @param[in] side  - Whether A is on the left ('L' or 'l') or right ('R', 'r').
    !> @param[in] uplo  - specifies  whether  the  upper ('U' or 'u')  or  lower ('L' or 'l') triangular  part  of  the  hermitian matrix is used
    !> @param[in] m     - the number of rows of C
    !> @param[in] n     - the number of columns of C
    !> @param[in] alpha - the \(\alpha\) scaling factor
    !> @param[in] da    - Host C-pointer to device memory to matrix A
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] db    - Host C-pointer to device memory to matrix B
    !> @param[in] ldb   - leading dimension of the matrix B
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] dc    - Host C-pointer to device memory to matrix C
    !> @param[in] ldc   - leading dimension of the matrix C
    !> @param[in] world - the device-host handler
    subroutine zhemm_gpu(side, uplo, m, n, alpha, da, lda, db, ldb, beta, dc, ldc, world)
        character, intent(in)               :: side
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r64), intent(in)            :: alpha
        type(c_ptr),  value                 :: da
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: db
        integer(i32), intent(in)            :: ldb
        complex(r64), intent(in)            :: beta
        type(c_ptr),  value                 :: dc
        integer(i32), intent(in)            :: ldc
        type(device_world_t), intent(inout) :: world

        integer(i32) :: ka
        integer(i32) :: magma_side, magma_uplo
        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)

#if defined(NVIDIAGPU) || defined(AMDGPU)

        magma_side = merge(MagmaLeft, MagmaRight, side == 'L' .or. side == 'l')
        magma_uplo = merge(MagmaLower, MagmaUpper, uplo == 'L' .or. uplo == 'l')
        
        call magma_zhemm(magma_side, magma_uplo, m, n, alpha, da, lda, db, ldb, beta, dc, ldc, world%get_queue())

#endif

#if defined(INTELGPU)

        ka = merge(m,n,side == 'L' .or. side == 'l')

        call c_f_pointer(da, A, [lda,ka])
        call c_f_pointer(db, B, [ldb,n])
        call c_f_pointer(dc, C, [ldc,n])

        !$omp dispatch is_device_ptr(a,b,c)
        call zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)

        nullify(A,B,C)
#endif

    end subroutine zhemm_gpu

    !> Computes dB = alpha*dA + beta*dB. Double precision
    !> @param[in] m     - The number of rows of the matrix A
    !> @param[in] n     - The number of columns of the matrix A
    !> @param[in] alpha - alpha scaling factor.
    !> @param[in] dA    - Host C-pointer to device A matrix.
    !> @param[in] lda   - The leading dimension of the array A.
    !> @param[in] beta  - beta scaling factor.
    !> @param[in,out] dB    - Host C-pointer to device B matrix.
    !> @param[in] ldb   - The leading dimension of the array B.
    !> @param[in,out] world - the device-host handler.
    subroutine zgeadd_gpu(m, n, alpha, dA, lda, beta, dB, ldb, world)
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r64), intent(in)            :: alpha
        type(c_ptr),  value                 :: dA
        integer(i32), intent(in)            :: lda
        complex(r64), intent(in)            :: beta
        type(c_ptr),  value                 :: dB
        integer(i32), intent(in)            :: ldb
        type(device_world_t), intent(inout) :: world

#if defined(INTELGPU)
        complex(r64), pointer :: A(:,:), B(:,:)
        call c_f_pointer(dA, A, [lda,n])
        call c_f_pointer(dB, B, [ldb,n])
        !$omp dispatch is_device_ptr(A,B)
        call mkl_zomatadd_batch_strided('c','n','n',m,n,alpha,A,lda,lda*m,beta,B,ldb,ldb*m,B,ldb,ldb*m,1)
        nullify(A,B)
#endif

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magmablas_zgeadd2(m, n, alpha, dA, lda, beta, dB, ldb, world%get_queue())
#endif

    end subroutine zgeadd_gpu

    !> Copies x contents to y.
    !> @param[in]       n           - Number of elements in vectors x and y. n >= 0.
    !> @param[in]       dx          - Host C-pointer to device memory to x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]       incx        - Stride between consecutive elements of dx.
    !> @param[in,out]   dy          - Host C-pointer to device memory to y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]       incy        - Stride between consecutive elements of dy.
    !> @param[in,out]   world       - the device-host handler
    subroutine zcopy_gpu(n, dx, incx, dy, incy, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zcopyvector_async(n, dx, incx, dy, incy, world%get_queue())
#endif

#if defined(INTELGPU)
        complex(r64), pointer :: x(:), y(:)
        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])
        !$omp dispatch is_device_ptr(x,y)
        call zcopy(n, x, incx, y, incy)
        nullify(x,y)
#endif

    end subroutine zcopy_gpu

    !> Sets an array to zero
    !> @param[in,out] dA    - Host pointer to device A array.
    !> @param[in]     size  - the size of the array
    !> @param[in,out] world - the device-host handler.
    subroutine zsetzero_gpu(dA, size, world)
        type(c_ptr),  value                 :: dA
        integer(c_size_t), intent(in)       :: size
        type(device_world_t), intent(inout) :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        integer(i32) :: ierr
        ierr =  magma_memset_async(dA, 0_c_int, size*bytes_double_complex, world%get_queue())
#endif

#if defined(INTELGPU)
        complex(r64), pointer   :: A(:)
        complex(r64), parameter :: zzero = cmplx(0.0_r64,0.0_r64,kind=r64) 
        integer(c_size_t) :: i
        call c_f_pointer(dA, A, [size])
        !$omp target has_device_addr(A)
        !$omp teams distribute parallel do simd
        do i = 1, size
            A(i) = zzero
        end do
        !$omp end teams distribute parallel do simd
        !$omp end target
        nullify(A) 
#endif

    end subroutine zsetzero_gpu 

end module device_linalg_common_interface
