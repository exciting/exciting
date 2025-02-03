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
    use omp_lib
    use iso_c_binding
    use magma2
#endif
#if defined(INTELGPU)
    use omp_lib
    use iso_c_binding
    use onemkl_lapack_omp_offload_lp64
    use onemkl_blas_omp_offload_lp64
#endif
    use iso_fortran_env,  only: i32=>int32, r32=>real32, r64=>real64
    use m_device_world_t, only: device_world_t

    implicit none

    private 
    public :: cgemm_gpu, cdotc_gpu, cdotu_gpu, cgetrf_gpu, cgetri_gpu, get_cgetri_nb_gpu, &
              zgemm_gpu, zdotc_gpu, zdotu_gpu, zgetrf_gpu, zgetri_gpu, get_zgetri_nb_gpu, &
              caxpy_gpu, zaxpy_gpu
              

contains

    !!!!!!!!!!!!!!!  SINGLE PRECISION !!!!!!!!!!!!!!

    !> Complex single precision LU decomposition.
    !> @param[in] m - The number of rows of the matrix A
    !> @param[in] n - The number of columns of the matrix A
    !> @param[in,out] dA - C pointer to A matrix (entry) in the device. On exit, is a C pointer to the factors L and U from the factorization.
    !> @param[in] lda - The leading dimension of the array A.
    !> @param[in,out] ipiv - The pivot indices; for 1 <= i <= min(m,n), row i of the matrix was interchanged with row ipiv(i). Host array.
    !> @param[out] info - execution information. 0 success; < 0:  if INFO = -i, the i-th argument had an illegal value;  > 0:  if INFO = i, U(i,i) is exactly zero.
    !> @param[in,out] world - the device-host handler
    subroutine cgetrf_gpu(m, n, dA, lda, ipiv, info, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        type(C_ptr),  value                             :: dA
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
        !$omp target variant dispatch use_device_ptr(A)
        call cgetrf(m, n, A, lda, ipiv, info)
        !$omp end target variant dispatch
        !$omp end target data

        nullify(A)
#endif
    end subroutine cgetrf_gpu

    !> Computes the inverse of a single precision complex matrix using its LU factorization.
    !> @param[in] n - The order of the matrix A
    !> @param[in,out] dA - On entry, the C-pointer to factors L and U from the factorization of matrix A. On exit, if successful, contains the address to the inverse of A. Device C-pointer.
    !> @param[in] ipiv - pivot indices. Host array.
    !> @param[in,out] dwork - Device pointer to work space.
    !> @param[in,out] lwork - size of the workspace
    !> @param[in,out] world - the device-host handler
    subroutine cgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(C_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(C_ptr),  value                             :: dwork
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
        !$omp target variant dispatch use_device_ptr(A)
        call cgetri(n, A, lda, ipiv, work, lwork, info)
        !$omp end target variant dispatch
        !$omp end target data

        nullify(A, work)
#endif

    end subroutine cgetri_gpu

    !> Provides the appropiate size for the workspace of the single precision complex matrix inversion using LU decomposition
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
    !> @param[in] da     - C-pointer to the A matrix. Device pointer
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] db     - C-pointer to the B matrix. Device pointer
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] dc - C-pointer to the C matrix. Device pointer
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in,out] world - device-host handler.
    subroutine cgemm_gpu(transa, transb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world)
        character, intent(in)         :: transa
        character, intent(in)         :: transb
        integer(i32), intent(in)      :: m
        integer(i32), intent(in)      :: n
        integer(i32), intent(in)      :: k
        complex(r32), intent(in)      :: alpha
        type(c_ptr),  value           :: da
        integer(i32), intent(in)      :: lda
        type(c_ptr),  value           :: db
        integer(i32), intent(in)      :: ldb
        complex(r32), intent(in)      :: beta
        type(c_ptr),  value           :: dc
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

        call magma_cgemm(opA, opB, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)
        integer(c_int) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(da, A, [lda,ka])
        call c_f_pointer(db, B, [ldb,kb])
        call c_f_pointer(dc, C, [ldc,n])

        !$omp target variant dispatch use_device_ptr(A,B,C)
        call cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        !$omp end target variant dispatch

        nullify(A, B, C)

#endif
    end subroutine cgemm_gpu

    !> Complex single precision dot product (unconjugated) of vectors x and y; \( x^T y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Device pointer
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Device pointer
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
        call c_f_pointer(dx, x, [n])
        call c_f_pointer(dy, y, [n])
        !$omp target data
        !$omp target variant dispatch
        cdotu_gpu = cdotu(n, x, incx, y, incy)
        !$omp end target variant dispatch
        !$omp end target data
        nullify(x, y)
#endif 
    end function cdotu_gpu

    !> Complex single precision dot product (conjugated) of vectors x and y; \( x^H y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Device pointer
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Device pointer
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
        call c_f_pointer(dx, x, [n])
        call c_f_pointer(dy, y, [n])
        
        !$omp target variant dispatch use_device_ptr(x,y)
        cdotc_gpu = cdotc(n, x, incx, y, incy)
        !$omp end target variant dispatch
        
        nullify(x, y)
#endif 
    end function cdotc_gpu

    !> Complex single precision constant times a vector plus a vector; \( y = \alpha x + y \). 
    !> @param[in]	n	- Number of elements in vectors x and y. n >= 0.
    !> @param[in]	alpha	- Scalar \( \alpha \)
    !> @param[in]	dx	- Device pointer to x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]	incx	- Stride between consecutive elements of dx. incx != 0.
    !> @param[in,out]	dy	- Device pointer to y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]	incy	- Stride between consecutive elements of dy. incy != 0.
    !> @param[in,out]   world	- the device-host handler
    subroutine caxpy_gpu(n, alpha, dx, incx, dy, incy, world)
        integer(i32), intent(in)                        :: n
        complex(r32), intent(in)                        :: alpha
        type(C_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(C_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_caxpy(n, alpha, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r32), pointer :: x(:), y(:)

        call c_f_pointer(dx, x, [1 + (n-1)*incx])
        call c_f_pointer(dy, y, [1 + (n-1)*incy])

        !$omp target variant dispatch use_device_ptr(x,y)
        call caxpy(n, alpha, x, incx, y, incy)
        !$omp end target variant dispatch

        nullify(x,y)

#endif 

    end subroutine caxpy_gpu

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
        type(C_ptr),  value                             :: dA
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
        !$omp target variant dispatch use_device_ptr(A)
        call zgetrf(m, n, A, lda, ipiv, info)
        !$omp end target variant dispatch
        !$omp end target data

        nullify(A)
#endif
    end subroutine zgetrf_gpu

    !> Computes the inverse of a double precision complex matrix using its LU factorization.
    !> @param[in] n - The order of the matrix A
    !> @param[in,out] dA - On entry, the C-pointer to factors L and U from the factorization of matrix A. On exit, if successful, contains the address to the inverse of A. Device C-pointer.
    !> @param[in] ipiv - pivot indices. Host array.
    !> @param[in,out] dwork - Device pointer to work space.
    !> @param[in,out] lwork - size of the workspace
    !> @param[in,out] world - the device-host handler
    subroutine zgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(C_ptr),  value                             :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(C_ptr),  value                             :: dwork
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
        !$omp target variant dispatch use_device_ptr(A)
        call zgetri(n, A, lda, ipiv, work, lwork, info)
        !$omp end target variant dispatch
        !$omp end target data

        nullify(A, work)
#endif

    end subroutine zgetri_gpu

    !> Provides the appropiate size for the workspace of the double precision complex matrix inversion using LU decomposition
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
    !> @param[in] da     - C-pointer to the A matrix. Device pointer
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] db     - C-pointer to the B matrix. Device pointer
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] dc - C-pointer to the C matrix. Device pointer
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in,out] world - device-host handler.
    subroutine zgemm_gpu(transa, transb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world)
        character, intent(in)            :: transa
        character, intent(in)            :: transb
        integer(i32), intent(in)         :: m
        integer(i32), intent(in)         :: n
        integer(i32), intent(in)         :: k
        complex(r64), intent(in)         :: alpha
        type(c_ptr),  value              :: da
        integer(i32), intent(in)         :: lda
        type(c_ptr),  value              :: db
        integer(i32), intent(in)         :: ldb
        complex(r64), intent(in)         :: beta
        type(c_ptr),  value              :: dc
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

        call magma_zgemm(opA, opB, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)
        integer(c_int) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(da, A, [lda,ka])
        call c_f_pointer(db, B, [ldb,kb])
        call c_f_pointer(dc, C, [ldc,n])

        !$omp target variant dispatch use_device_ptr(a,b,c) 
        call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        !$omp end target variant dispatch

        nullify(A, B, C)

#endif
    end subroutine zgemm_gpu

    !> Complex double precision dot product (unconjugated) of vectors x and y; \( x^T y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Device pointer
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Device pointer
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
        call c_f_pointer(dx, x, [n])
        call c_f_pointer(dy, y, [n])
        !$omp target data
        !$omp target variant dispatch
        zdotu_gpu = zdotu(n, x, incx, y, incy)
        !$omp end target variant dispatch
        !$omp end target data
        nullify(x, y)
#endif 
    end function zdotu_gpu

    !> Complex double precision dot product (conjugated) of vectors x and y; \( x^H y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] dx - C-pointer to the x vector. Device pointer
    !> @param[in] incx - Stride between consecutive elements of dx
    !> @param[in] dy - C-pointer to the y vector. Device pointer
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
        call c_f_pointer(dx, x, [n])
        call c_f_pointer(dy, y, [n])
        
        !$omp target variant dispatch use_device_ptr(x,y)
        zdotc_gpu = zdotc(n, x, incx, y, incy)
        !$omp end target variant dispatch
        
        nullify(x, y)
#endif 
    end function zdotc_gpu

    !> Complex double precision constant times a vector plus a vector; \( y = \alpha x + y \). 
    !> @param[in]	n	- Number of elements in vectors x and y. n >= 0.
    !> @param[in]	alpha	- Scalar \( \alpha \)
    !> @param[in]	dx	- Device pointer to x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]	incx	- Stride between consecutive elements of dx. incx != 0.
    !> @param[in,out]	dy	- Device pointer to y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]	incy	- Stride between consecutive elements of dy. incy != 0.
    !> @param[in,out]   world	- the device-host handler
    subroutine zaxpy_gpu(n, alpha, dx, incx, dy, incy, world)
        integer(i32), intent(in)                        :: n
        complex(r64), intent(in)                        :: alpha
        type(C_ptr),  value                             :: dx
        integer(i32), intent(in)                        :: incx
        type(C_ptr),  value                             :: dy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zaxpy(n, alpha, dx, incx, dy, incy, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r64), pointer :: x(:), y(:)

        call c_f_pointer(dx, x, [incx*n])
        call c_f_pointer(dy, y, [incy*n])

        !$omp target variant dispatch use_device_ptr(x,y)
        call zaxpy(n, alpha, x, incx, y, incy)
        !$omp end target variant dispatch

        nullify(x,y)

#endif 

    end subroutine zaxpy_gpu

end module device_linalg_common_interface
