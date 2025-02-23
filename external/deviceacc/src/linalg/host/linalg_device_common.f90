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
!> This is the CPU backend for compilations not supporting offload

!> Module containing unified calls to linear algebra device accelerated routines
!> CPU backend; i.e. no device accelerated routines are called
module device_linalg_common_interface
        
    use iso_c_binding,    only: c_ptr, c_int, c_size_t, c_f_pointer
    use iso_fortran_env,  only: i32=>int32, r32=>real32, r64=>real64
    use m_memory_device,  only: bytes_single_complex, bytes_double_complex, generate_batched_array
    use m_device_world_t, only: device_world_t

    implicit none

    private 
    public :: cgemm_gpu, cdotc_gpu, cdotu_gpu, cgetrf_gpu, cgetri_gpu, get_cgetri_nb_gpu, &
              zgemm_gpu, zdotc_gpu, zdotu_gpu, zgetrf_gpu, zgetri_gpu, get_zgetri_nb_gpu, &
              caxpy_gpu, zaxpy_gpu, cgemm_batched_gpu, zgemm_batched_gpu, &
              cgerc_gpu, cgeru_gpu, zgerc_gpu, zgeru_gpu, &
              chemv_gpu, zhemv_gpu, chemm_gpu, zhemm_gpu, &
              cgeadd_gpu, zgeadd_gpu, &
              ccopy_gpu, zcopy_gpu, &
              csetzero_gpu, zsetzero_gpu  


   interface
       subroutine memset(ptr, value, size) bind(c, name="memset")
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: ptr
            integer(c_int), value :: value
            integer(c_size_t), value :: size
        end subroutine memset
   end interface

contains

!!!!!!!!!!!!!!!  SINGLE PRECISION !!!!!!!!!!!!!!

    !> Complex single precision LU decomposition.
    !> @param[in] m - The number of rows of the matrix A
    !> @param[in] n - The number of columns of the matrix A
    !> @param[in,out] hA - C pointer to A matrix (entry) in the host. On exit, is a C pointer to the factors L and U from the factorization.
    !> @param[in] lda - The leading dimension of the array A.
    !> @param[in,out] ipiv - The pivot indices; for 1 <= i <= min(m,n), row i of the matrix was interchanged with row ipiv(i). Host array.
    !> @param[out] info - execution information. 0 success; < 0:  if INFO = -i, the i-th argument had an illegal value;  > 0:  if INFO = i, U(i,i) is exactly zero.
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine cgetrf_gpu(m, n, hA, lda, ipiv, info, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(device_world_t), intent(inout)             :: world

        complex(r32), pointer :: A(:,:)

        call c_f_pointer(hA, A, [lda,n])
        call cgetrf(m, n, A, lda, ipiv, info)

        nullify(A)

    end subroutine cgetrf_gpu

    !> Computes the inverse of a single precision complex matrix using its LU factorization.
    !> @param[in] n - The order of the matrix A
    !> @param[in,out] hA - On entry, the C-pointer to factors L and U from the factorization of matrix A. On exit, if successful,
    !>                contains the address to the inverse of A. Host C-pointer to host memory.
    !> @param[in] ipiv - pivot indices. Host array.
    !> @param[in,out] dwork - Host pointer to work space.
    !> @param[in,out] lwork - size of the workspace
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine cgetri_gpu(n, hA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(c_ptr),  value                             :: dwork
        integer(i32), intent(in)                        :: lwork
        type(device_world_t), intent(inout)             :: world

        complex(r32), pointer :: A(:,:), work(:)

        call c_f_pointer(hA, A, [lda,n])
        call c_f_pointer(dwork, work, [lwork])
        call cgetri(n, A, lda, ipiv, work, lwork, info)
        nullify(A, work)

    end subroutine cgetri_gpu

    !> Provides the appropriate block size for the workspace of the single precision complex matrix inversion using LU decomposition
    !> The optimal worksize is then obtained by multiplying n by the result of this function
    !> @param[in] n - the order of the A matrix
    pure function get_cgetri_nb_gpu(n) result(nblock)
        integer(i32), intent(in) :: n
        integer :: nblock
        nblock = 64_i32
    end function get_cgetri_nb_gpu

    !> Complex single precision matrix-matrix product (C = alpha * op(A) * op(B) + beta * C)
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      -  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] hA     - C-pointer to the A matrix. Host C-pointer to host memory.
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] hB     - C-pointer to the B matrix. Host C-pointer to host memory.
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] hC - C-pointer to the C matrix. Host C-pointer to host memory.
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in,out] world - device-host handler. CPU backend.
    subroutine cgemm_gpu(transa, transb, m, n, k, alpha, hA, lda, hB, ldb, beta, hC, ldc, world)
        character, intent(in)               :: transa
        character, intent(in)               :: transb
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        integer(i32), intent(in)            :: k
        complex(r32), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: hB
        integer(i32), intent(in)            :: ldb
        complex(r32), intent(in)            :: beta
        type(c_ptr),  value                 :: hC
        integer(i32), intent(in)            :: ldc
        type(device_world_t), intent(inout) :: world

        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)
        integer(i32) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(hA, A, [lda,ka])
        call c_f_pointer(hB, B, [ldb,kb])
        call c_f_pointer(hC, C, [ldc,n])

        call cgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

        nullify(A, B, C)

    end subroutine cgemm_gpu

    !> Complex single precision matrix-matrix batched product (C[ibatch] = alpha * op(A[ibatch]) * op(B[ibatch]) + beta * C[ibatch])
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      -  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] hA     - C-pointer to the A matrices. CPU-pointer to host memory
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] hB     - C-pointer to the B matrices. CPU-pointer to host memory
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] hC - C-pointer to the C matrices. CPU-pointer to host memory
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in] batchcount - number of matrices in the batch
    !> @param[in,out] world - device-host handler. CPU backend.
    !> @param[in,optional] stridea - If present it indicates a different stride for A matrix than the matrix size.
    !> @param[in,optional] strideb - If present it indicates a different stride for B matrix than the matrix size.
    !> @param[in,optional] stridec - If present it indicates a different stride for C matrix than the matrix size.
    subroutine cgemm_batched_gpu(transa, transb, m, n, k, alpha, hA, lda, hB, ldb, beta, hC, ldc, batchcount, world, stridea, &
                                 strideb, stridec)
        character, intent(in)            :: transa
        character, intent(in)            :: transb
        integer(i32), intent(in)         :: m
        integer(i32), intent(in)         :: n
        integer(i32), intent(in)         :: k
        complex(r32), intent(in)         :: alpha
        type(c_ptr),  value              :: hA
        integer(i32), intent(in)         :: lda
        type(c_ptr),  value              :: hB
        integer(i32), intent(in)         :: ldb
        complex(r32), intent(in)         :: beta
        type(c_ptr),  value              :: hC
        integer(i32), intent(in)         :: ldc
        integer(i32), intent(in)         :: batchcount
        type(device_world_t), intent(in) :: world
        integer(i32), intent(in), optional :: stridea
        integer(i32), intent(in), optional :: strideb
        integer(i32), intent(in), optional :: stridec

        integer(c_int) :: ka, kb
        integer(c_int) :: ibatch

        type(c_ptr), target :: hA_array(batchcount), hB_array(batchcount), hC_array(batchcount)
        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)
        integer(i32) :: stridea_local, strideb_local, stridec_local

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

#if defined(_INTEL_COMPILER)

        call c_f_pointer(hA, A, [lda,ka])
        call c_f_pointer(hB, B, [ldb,kb])
        call c_f_pointer(hC, C, [ldc,n])
        call cgemm_batched_strided(transa, transb, m, n, k, alpha, A, lda, stridea_local, B, ldb, &
                                   strideb_local, beta, C, ldc, stridec_local, batchcount)
        nullify(A, B, C)

#else

        call generate_batched_array(hA, int(batchcount,kind=c_size_t), int(stridea_local,kind=c_size_t)*bytes_single_complex, hA_array)
        call generate_batched_array(hB, int(batchcount,kind=c_size_t), int(strideb_local,kind=c_size_t)*bytes_single_complex, hB_array)
        call generate_batched_array(hC, int(batchcount,kind=c_size_t), int(stridec_local,kind=c_size_t)*bytes_single_complex, hC_array)

        !$omp parallel do &
        !$omp default(none) shared(lda,ka,ldb,kb,ldc,n,transa,transb,m,k,alpha,beta,batchcount,hA_array,hB_array,hC_array) &
        !$omp private(ibatch,A,B,C)
        do ibatch = 1, batchcount
            call c_f_pointer(hA_array(ibatch), A, [lda,ka])
            call c_f_pointer(hB_array(ibatch), B, [ldb,kb])
            call c_f_pointer(hC_array(ibatch), C, [ldc,n])
            call cgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            nullify(A, B, C)
        end do
        !$omp end parallel do

#endif

    end subroutine cgemm_batched_gpu


    !> Complex single precision dot product (unconjugated) of vectors x and y; \( x^T y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] hx - C-pointer to the x vector. Host pointer to host memory
    !> @param[in] incx - Stride between consecutive elements of hx
    !> @param[in] hy - C-pointer to the y vector. Host pointer to host memory
    !> @param[in] incy - Stride between consecutive elements of hy
    !> @param[in,out] world - device-host handler.
    complex(r32) function cdotu_gpu(n, hx, incx, hy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: hx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: hy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world

        complex(r32), external :: cdotu
 
        complex(r32), pointer :: x(:), y(:)

        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        cdotu_gpu = cdotu(n, x, incx, y, incy)
        nullify(x, y)
    end function cdotu_gpu

    !> Complex single precision dot product (conjugated) of vectors x and y; \( x^H y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] hx - C-pointer to the x vector. Host pointer to host memory
    !> @param[in] incx - Stride between consecutive elements of hx
    !> @param[in] hy - C-pointer to the y vector. Host pointer to host memory
    !> @param[in] incy - Stride between consecutive elements of hy
    !> @param[in,out] world - device-host handler.
    complex(r32) function cdotc_gpu(n, hx, incx, hy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: hx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: hy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world

        complex(r32), external :: cdotc
        
        complex(r32), pointer :: x(:), y(:)
        
        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        cdotc_gpu = cdotc(n, x, incx, y, incy)
        nullify(x, y)

    end function cdotc_gpu

    !> Complex single precision constant times a vector plus a vector; \( y = \alpha x + y \). 
    !> @param[in]    	n	        - Number of elements in vectors x and y. n >= 0.
    !> @param[in]	    alpha	    - Scalar \( \alpha \)
    !> @param[in]	    hx	        - Host pointer to x to host memory. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]	    incx	    - Stride between consecutive elements of hx. incx != 0.
    !> @param[in,out]	hy	        - Host pointer to y to host memory. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]	    incy	    - Stride between consecutive elements of hy. incy != 0.
    !> @param[in,out]   world	    - the device-host handler (CPU-backend)
    subroutine caxpy_gpu(n, alpha, hx, incx, hy, incy, world)
        integer(i32), intent(in)                        :: n
        complex(r32), intent(in)                        :: alpha
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

        complex(r32), pointer :: x(:), y(:)

        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call caxpy(n, alpha, x, incx, y, incy)
        nullify(x,y)

    end subroutine caxpy_gpu


    !> Perform rank-1 update, \( A = \alpha x y^H + A \) for complex single precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] hx -    Host pointer to vector x to host memory. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] hy - Host pointer to vector y to host memory. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y to host memory. Must not be zero.
    !> @param[in,out] hA - Host pointer to A. On entry, the m-by-n matrix A. On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldhA >= max(1, m))
    !> @param[in,out]   world - the device-host handler (CPU-backend)
    subroutine cgerc_gpu(m, n, alpha, hx, incx, hy, incy, hA, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r32), intent(in)                        :: alpha
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world
        
        complex(r32), pointer :: x(:), y(:), A(:,:)

        call c_f_pointer(hx, x, [1 + (m-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call c_f_pointer(hA, A, [lda,n])
        call cgerc(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)

    end subroutine cgerc_gpu

    !> Perform rank-1 update, \( A = \alpha x y^T + A \) for complex single precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] hx -    Host pointer to vector x to host memory. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] hy - Host pointer to vector y to host memory. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y. Must not be zero.
    !> @param[in,out] hA - Host pointer to A to host memory. On entry, the m-by-n matrix A. 
    !>                On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldhA >= max(1, m))
    !> @param[in,out]   world - the device-host handler (CPU-backend)
    subroutine cgeru_gpu(m, n, alpha, hx, incx, hy, incy, hA, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r32), intent(in)                        :: alpha
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world
        
        complex(r32), pointer :: x(:), y(:), A(:,:)

        call c_f_pointer(hx, x, [1 + (m-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call c_f_pointer(hA, A, [lda,n])
        call cgeru(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)

    end subroutine cgeru_gpu


    !> Perform Hermitian matrix-vector product in single precision.
    !> \( y = \alpha A x + \beta y \) (side == MagmaLeft), or
    !> \( y = \alpha A v + \beta y \) (side == MagmaRight),
    !> where \( A \) is Hermitian.
    !>
    !> @param[in] uplo  - specifies  whether  the  upper ('U' or 'u')  or  lower ('L' or 'l') triangular  part  of  the  hermitian matrix is used
    !> @param[in] n     - the order of the A matrix
    !> @param[in] alpha - the \(\alpha\) scaling factor
    !> @param[in] hA    - Host pointer to matrix A (host memory)
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] hx    - Host pointer to vector x (host memory)
    !> @param[in] incx  - Specifies the increment for the elements of x.
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] hy    - Host pointer to vector y (host memory)
    !> @param[in] incy  - Specifies the increment for the elements of y.
    !> @param[in] world - the device-host handler
    subroutine chemv_gpu(uplo, n, alpha, hA, lda, hx, incx, beta, hy, incy, world)
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: n
        complex(r32), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: hx
        integer(i32), intent(in)            :: incx
        complex(r32), intent(in)            :: beta
        type(c_ptr),  value                 :: hy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world

        complex(r32), pointer :: A(:,:), x(:), y(:)

        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call c_f_pointer(hA, A, [lda,n])

        call chemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)

        nullify(A,x,y)

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
    !> @param[in] hA    - Host pointer to matrix A (host memory)
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] hB    - Host pointer to matrix B (host memory)
    !> @param[in] ldb   - leading dimension of the matrix B
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] hC    - Host pointer to matrix C (host memory)
    !> @param[in] ldc   - leading dimension of the matrix C
    !> @param[in] world - the device-host handler. CPU-only backend
    subroutine chemm_gpu(side, uplo, m, n, alpha, hA, lda, hB, ldb, beta, hC, ldc, world)
        character, intent(in)               :: side
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r32), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: hB
        integer(i32), intent(in)            :: ldb
        complex(r32), intent(in)            :: beta
        type(c_ptr),  value                 :: hC
        integer(i32), intent(in)            :: ldc
        type(device_world_t), intent(inout) :: world

        integer(i32) :: ka
        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)

        ka = merge(m,n,side == 'L' .or. side == 'l')

        call c_f_pointer(hA, A, [lda,ka])
        call c_f_pointer(hB, B, [ldb,n])
        call c_f_pointer(hC, C, [ldc,n])
        call chemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        nullify(A,B,C)

    end subroutine chemm_gpu

    !> Computes hB = alpha*hA + beta*hB. Single precision 
    !> @param[in] m     - The number of rows of the matrix A
    !> @param[in] n     - The number of columns of the matrix A
    !> @param[in] alpha - alpha scaling factor.
    !> @param[in] hA    - Host c-pointer to host A matrix (entry).
    !> @param[in] lda   - The leading dimension of the array A.
    !> @param[in] beta  - beta scaling factor.
    !> @param[in,out] hB    - Host c-pointer to host B matrix. On exit it have the result.
    !> @param[in] ldb   - The leading dimension of the array B.
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine cgeadd_gpu(m, n, alpha, hA, lda, beta, hB, ldb, world)
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r32), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        complex(r32), intent(in)            :: beta
        type(c_ptr),  value                 :: hB
        integer(i32), intent(in)            :: ldb
        type(device_world_t), intent(inout) :: world

        complex(r32), pointer :: A(:,:), B(:,:)
        integer(i32) :: i, j

        call c_f_pointer(hA, A, [lda,n])
        call c_f_pointer(hB, B, [ldb,n])

#if defined(_INTEL_COMPILER) 
        call mkl_comatadd('c','n','n',m,n,A,lda,B,ldb,B,ldb)
#else
        !$omp parallel do collapse(2) default(shared) private(i,j)
        do i = 1, n
            do j = 1, m
                B(j,i) = beta * B(j,i) + alpha * A(j,i)
            end do
        end do
        !$omp end parallel do
#endif

        nullify(A,B)

    end subroutine cgeadd_gpu

    !> Copies x contents to y.
    !> @param[in]       n           - Number of elements in vectors x and y. n >= 0.
    !> @param[in]       hx          - Host pointer to host x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]       incx        - Stride between consecutive elements of hx.
    !> @param[in,out]   hy          - Host pointer to host y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]       incy        - Stride between consecutive elements of hy.
    !> @param[in,out]   world       - the device-host handler (CPU-backend)
    subroutine ccopy_gpu(n, hx, incx, hy, incy, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

        complex(r32), pointer :: x(:), y(:)

        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call ccopy(n, x, incx, y, incy)
        nullify(x,y)

    end subroutine ccopy_gpu

    !> Sets an array to zero
    !> @param[in,out] hA    - Host pointer to host A array.
    !> @param[in,out] size  - the size of the array
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine csetzero_gpu(hA, size, world)
        type(c_ptr),  value                 :: hA
        integer(c_size_t), intent(in)       :: size
        type(device_world_t), intent(inout) :: world

        call memset(hA, 0_c_int, size*bytes_single_complex)

    end subroutine csetzero_gpu

    ! DOUBLE PRECISION

    !> Complex double precision LU decomposition.
    !> @param[in] m - The number of rows of the matrix A
    !> @param[in] n - The number of columns of the matrix A
    !> @param[in,out] hA - C pointer to A matrix (entry) in the host. On exit, is a C pointer to the factors L and U from the factorization.
    !> @param[in] lda - The leading dimension of the array A.
    !> @param[in,out] ipiv - The pivot indices; for 1 <= i <= min(m,n), row i of the matrix was interchanged with row ipiv(i). Host array.
    !> @param[out] info - execution information. 0 success; < 0:  if INFO = -i, the i-th argument had an illegal value;  > 0:  if INFO = i, U(i,i) is exactly zero.
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine zgetrf_gpu(m, n, hA, lda, ipiv, info, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(device_world_t), intent(inout)             :: world

        complex(r64), pointer :: A(:,:)
        call c_f_pointer(hA, A, [lda,n])
        call zgetrf(m, n, A, lda, ipiv, info)
        nullify(A)

    end subroutine zgetrf_gpu

    !> Computes the inverse of a double precision complex matrix using its LU factorization.
    !> @param[in] n - The order of the matrix A
    !> @param[in,out] hA - On entry, the C-pointer to factors L and U from the factorization of matrix A. On exit, if successful,
    !>                contains the address to the inverse of A. Host C-pointer (host memory).
    !> @param[in] ipiv - pivot indices. Host array.
    !> @param[in,out] dwork - Host pointer to work space.
    !> @param[in,out] lwork - size of the workspace
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine zgetri_gpu(n, hA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(c_ptr),  value                             :: dwork
        integer(i32), intent(in)                        :: lwork
        type(device_world_t), intent(inout)             :: world

        complex(r64), pointer :: A(:,:), work(:)

        call c_f_pointer(hA, A, [lda,n])
        call c_f_pointer(dwork, work, [lwork])
        call zgetri(n, A, lda, ipiv, work, lwork, info)
        nullify(A, work)

    end subroutine zgetri_gpu

    !> Provides the appropriate block size for the workspace of the double precision complex matrix inversion using LU decomposition
    !> The optimal worksize is then obtained by multiplying n by the result of this function
    !> @param[in] n - the order of the A matrix
    pure function get_zgetri_nb_gpu(n) result(nblock)
        integer(i32), intent(in) :: n
        integer :: nblock
        nblock = 64_i32
    end function get_zgetri_nb_gpu

    !> Complex double precision matrix-matrix product (C = alpha * op(A) * op(B) + beta * C)
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      - the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] hA     - C-pointer to the A matrix. Host pointer to host memory
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] hB     - C-pointer to the B matrix. Host pointer to host memory
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] hC - C-pointer to the C matrix. Host pointer to host memory
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in,out] world - device-host handler. CPU backend.
    subroutine zgemm_gpu(transa, transb, m, n, k, alpha, hA, lda, hB, ldb, beta, hC, ldc, world)
        character, intent(in)               :: transa
        character, intent(in)               :: transb
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        integer(i32), intent(in)            :: k
        complex(r64), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: hB
        integer(i32), intent(in)            :: ldb
        complex(r64), intent(in)            :: beta
        type(c_ptr),  value                 :: hC
        integer(i32), intent(in)            :: ldc
        type(device_world_t), intent(inout) :: world

        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)
        integer(i32) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(hA, A, [lda,ka])
        call c_f_pointer(hB, B, [ldb,kb])
        call c_f_pointer(hC, C, [ldc,n])

        call zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

        nullify(A, B, C)

    end subroutine zgemm_gpu

    !> Complex single precision matrix-matrix batched product (C[ibatch] = alpha * op(A[ibatch]) * op(B[ibatch]) + beta * C[ibatch])
    !> @param[in] transa - specifies the form of op( A ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] transb - specifies the form of op( B ) to be used in the matrix multiplication as follows: 'n'/'N' nothing; 'T'/'t' transpose; 'C'/'c' adjoint.
    !> @param[in] m      - specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.
    !> @param[in] n      - specifies  the number  of columns  of the  matrix op( B )  and of the  matrix  C.
    !> @param[in] k      -  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
    !> @param[in] alpha  - specifies the scalar alpha multiplying the op(A) * op(B)
    !> @param[in] hA     - C-pointer to the A matrices. CPU-pointer to host memory
    !> @param[in] lda    - leading dimension of A matrix
    !> @param[in] hB     - C-pointer to the B matrices. CPU-pointer to host memory
    !> @param[in] ldb    - leading dimension of B matrix
    !> @param[in] beta   - specifies the scalar  beta so that a fraction of the old C is added to the matrix product alpha * op(A) * op(B)
    !> @param[in,out] hC - C-pointer to the C matrices. CPU-pointer to host memory
    !> @param[in] ldc    - leading dimension of C matrix
    !> @param[in] batchcount - number of matrices in the batch
    !> @param[in,out] world - device-host handler. CPU backend.
    !> @param[in,optional] stridea - If present it indicates a different stride for A matrix than the matrix size.
    !> @param[in,optional] strideb - If present it indicates a different stride for B matrix than the matrix size.
    !> @param[in,optional] stridec - If present it indicates a different stride for C matrix than the matrix size.
    subroutine zgemm_batched_gpu(transa, transb, m, n, k, alpha, hA, lda, hB, ldb, beta, hC, ldc, batchcount, world, stridea, &
                                 strideb, stridec)
        character, intent(in)            :: transa
        character, intent(in)            :: transb
        integer(i32), intent(in)         :: m
        integer(i32), intent(in)         :: n
        integer(i32), intent(in)         :: k
        complex(r64), intent(in)         :: alpha
        type(c_ptr),  value              :: hA
        integer(i32), intent(in)         :: lda
        type(c_ptr),  value              :: hB
        integer(i32), intent(in)         :: ldb
        complex(r64), intent(in)         :: beta
        type(c_ptr),  value              :: hC
        integer(i32), intent(in)         :: ldc
        integer(i32), intent(in)         :: batchcount
        type(device_world_t), intent(in) :: world
        integer(i32), intent(in), optional :: stridea
        integer(i32), intent(in), optional :: strideb
        integer(i32), intent(in), optional :: stridec

        integer(c_int) :: ka, kb
        integer(c_int) :: ibatch

        type(c_ptr), target :: hA_array(batchcount), hB_array(batchcount), hC_array(batchcount)
        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)
        integer(i32) :: stridea_local, strideb_local, stridec_local

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

#if defined(_INTEL_COMPILER)

        call c_f_pointer(hA, A, [lda,ka])
        call c_f_pointer(hB, B, [ldb,kb])
        call c_f_pointer(hC, C, [ldc,n])
        call zgemm_batched_strided(transa, transb, m, n, k, alpha, A, lda, stridea_local, B, ldb, &
                                   strideb_local, beta, C, ldc, stridec_local, batchcount)
        nullify(A, B, C)

#else

        call generate_batched_array(hA, int(batchcount,kind=c_size_t), int(stridea_local,kind=c_size_t)*bytes_double_complex, hA_array)
        call generate_batched_array(hB, int(batchcount,kind=c_size_t), int(strideb_local,kind=c_size_t)*bytes_double_complex, hB_array)
        call generate_batched_array(hC, int(batchcount,kind=c_size_t), int(stridec_local,kind=c_size_t)*bytes_double_complex, hC_array)

        !$omp parallel do &
        !$omp default(none) shared(lda,ka,ldb,kb,ldc,n,transa,transb,m,k,alpha,beta,batchcount,hA_array,hB_array,hC_array) &
        !$omp private(ibatch,A,B,C)
        do ibatch = 1, batchcount
            call c_f_pointer(hA_array(ibatch), A, [lda,ka])
            call c_f_pointer(hB_array(ibatch), B, [ldb,kb])
            call c_f_pointer(hC_array(ibatch), C, [ldc,n])
            call zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            nullify(A, B, C)
        end do
        !$omp end parallel do

#endif

    end subroutine zgemm_batched_gpu

    !> Complex double precision dot product (unconjugated) of vectors x and y; \( x^T y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] hx - C-pointer to the host x vector.
    !> @param[in] incx - Stride between consecutive elements of hx
    !> @param[in] hy - C-pointer to the host y vector. 
    !> @param[in] incy - Stride between consecutive elements of hy
    !> @param[in,out] world - device-host handler.
    complex(r64) function zdotu_gpu(n, hx, incx, hy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: hx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: hy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world
 
        complex(r64), external :: zdotu

        complex(r64), pointer :: x(:), y(:)
        
        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        zdotu_gpu = zdotu(n, x, incx, y, incy)
        nullify(x, y)

    end function zdotu_gpu

    !> Complex double precision dot product (conjugated) of vectors x and y; \( x^H y \).
    !> @param[in] n - number of elements in vector x and y
    !> @param[in] hx - C-pointer to the host x vector.
    !> @param[in] incx - Stride between consecutive elements of hx
    !> @param[in] hy - C-pointer to the host y vector.
    !> @param[in] incy - Stride between consecutive elements of hy
    !> @param[in,out] world - device-host handler.
    complex(r64) function zdotc_gpu(n, hx, incx, hy, incy, world)
        integer(i32), intent(in)            :: n
        type(c_ptr),  value                 :: hx
        integer(i32), intent(in)            :: incx
        type(c_ptr),  value                 :: hy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world
        
        complex(r64), external :: zdotc

        complex(r64), pointer :: x(:), y(:)
        
        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        zdotc_gpu = zdotc(n, x, incx, y, incy)
        nullify(x, y)

    end function zdotc_gpu

    !> Complex double precision constant times a vector plus a vector; \( y = \alpha x + y \). 
    !> @param[in]    	n	        - Number of elements in vectors x and y. n >= 0.
    !> @param[in]	    alpha	    - Scalar \( \alpha \)
    !> @param[in]	    hx	        - Host pointer to host x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]	    incx	    - Stride between consecutive elements of hx. incx != 0.
    !> @param[in,out]	hy	        - Host pointer to host y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]	    incy	    - Stride between consecutive elements of hy. incy != 0.
    !> @param[in,out]   world	    - the device-host handler (CPU-backend)
    subroutine zaxpy_gpu(n, alpha, hx, incx, hy, incy, world)
        integer(i32), intent(in)                        :: n
        complex(r64), intent(in)                        :: alpha
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world
        
        complex(r64), pointer :: x(:), y(:)

        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call zaxpy(n, alpha, x, incx, y, incy)
        nullify(x,y)

    end subroutine zaxpy_gpu

    !> Perform rank-1 update, \( A = \alpha x y^H + A \) for complex double precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] hx -    Host pointer to host vector x. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] hy - Host pointer to vector y. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y. Must not be zero.
    !> @param[in,out] hA - Host pointer to host A. On entry, the m-by-n matrix A. On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldhA >= max(1, m))
    !> @param[in,out]   world - the device-host handler (CPU-backend)
    subroutine zgerc_gpu(m, n, alpha, hx, incx, hy, incy, hA, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r64), intent(in)                        :: alpha
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world
        
        complex(r64), pointer :: x(:), y(:), A(:,:)

        call c_f_pointer(hx, x, [1 + (m-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call c_f_pointer(hA, A, [lda,n])
        call zgerc(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)

    end subroutine zgerc_gpu

    !> Perform rank-1 update, \( A = \alpha x y^T + A \) for complex double precision arrays
    !> @param[in] m - The number of rows of the matrix A and the length of vector x. (m >= 0)
    !> @param[in] n - The number of columns of the matrix A and the length of vector y. (n >= 0)
    !> @param[in] alpha - The scalar multiplier.
    !> @param[in] hx -    Host pointer to vector x. This vector must be of size at least (1 + (m - 1) * abs(incx)).
    !> @param[in] incx - The increment for the elements of x. Must not be zero.
    !> @param[in] hy - Host pointer to vector y. This vector must be of size at least (1 + (n - 1) * abs(incy)).
    !> @param[in] incy - The increment for the elements of y. Must not be zero.
    !> @param[in,out] hA - Host pointer to A. On entry, the m-by-n matrix A. On exit, the updated matrix A after the rank-1 update.
    !> @param[in] lda - The leading dimension of the array A. (ldhA >= max(1, m))
    !> @param[in,out]   world - the device-host handler (CPU-backend)
    subroutine zgeru_gpu(m, n, alpha, hx, incx, hy, incy, hA, lda, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        complex(r64), intent(in)                        :: alpha
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(c_ptr),  value                             :: hA
        integer(i32), intent(in)                        :: lda
        type(device_world_t), intent(inout)             :: world
        
        complex(r64), pointer :: x(:), y(:), A(:,:)

        call c_f_pointer(hx, x, [1 + (m-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call c_f_pointer(hA, A, [lda,n])
        call zgeru(m, n, alpha, x, incx, y, incy, A, lda)
        nullify(x,y,A)

    end subroutine zgeru_gpu

    !> Perform Hermitian matrix-vector product in double precision.
    !> \( y = \alpha A x + \beta y \) (side == MagmaLeft), or
    !> \( y = \alpha A v + \beta y \) (side == MagmaRight),
    !> where \( A \) is Hermitian.
    !>
    !> @param[in] uplo  - specifies  whether  the  upper ('U' or 'u')  or  lower ('L' or 'l') triangular  part  of  the  hermitian matrix is used
    !> @param[in] n     - the order of the A matrix
    !> @param[in] alpha - the \(\alpha\) scaling factor
    !> @param[in] hA    - Host pointer to matrix A
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] hx    - Host pointer to vector x
    !> @param[in] incx  - Specifies the increment for the elements of x.
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] hy    - Host pointer to vector y
    !> @param[in] incy  - Specifies the increment for the elements of y.
    !> @param[in] world - the device-host handler
    subroutine zhemv_gpu(uplo, n, alpha, hA, lda, hx, incx, beta, hy, incy, world)
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: n
        complex(r64), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: hx
        integer(i32), intent(in)            :: incx
        complex(r64), intent(in)            :: beta
        type(c_ptr),  value                 :: hy
        integer(i32), intent(in)            :: incy
        type(device_world_t), intent(inout) :: world

        complex(r64), pointer :: A(:,:), x(:), y(:)

        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call c_f_pointer(hA, A, [lda,n])

        call zhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)

        nullify(A,x,y)

    end subroutine zhemv_gpu

    !> Perform Hermitian matrix-matrix product in double precision.
    !> \( C = \alpha A B + \beta C \) (side == MagmaLeft), or
    !> \( C = \alpha B A + \beta C \) (side == MagmaRight),
    !> where \( A \) is Hermitian.
    !>
    !> @param[in] side  - Whether A is on the left ('L' or 'l') or right ('R', 'r').
    !> @param[in] uplo  - specifies  whether  the  upper ('U' or 'u')  or  lower ('L' or 'l') triangular  part  of  the  hermitian matrix is used
    !> @param[in] m     - the number of rows of C
    !> @param[in] n     - the number of columns of C
    !> @param[in] alpha - the \(\alpha\) scaling factor
    !> @param[in] hA    - Host pointer to matrix A
    !> @param[in] lda   - leading dimension of the matrix A
    !> @param[in] hB    - Host pointer to matrix B
    !> @param[in] ldb   - leading dimension of the matrix B
    !> @param[in] beta  - the \(\beta\) scaling factor
    !> @param[in] hC    - Host pointer to matrix C
    !> @param[in] ldc   - leading dimension of the matrix C
    !> @param[in] world - the device-host handler. CPU-only backend
    subroutine zhemm_gpu(side, uplo, m, n, alpha, hA, lda, hB, ldb, beta, hC, ldc, world)
        character, intent(in)               :: side
        character, intent(in)               :: uplo
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r64), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        type(c_ptr),  value                 :: hB
        integer(i32), intent(in)            :: ldb
        complex(r64), intent(in)            :: beta
        type(c_ptr),  value                 :: hC
        integer(i32), intent(in)            :: ldc
        type(device_world_t), intent(inout) :: world

        integer(i32) :: ka
        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)

        ka = merge(m,n,side == 'L' .or. side == 'l')

        call c_f_pointer(hA, A, [lda,ka])
        call c_f_pointer(hB, B, [ldb,n])
        call c_f_pointer(hC, C, [ldc,n])
        call zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        nullify(A,B,C)

    end subroutine zhemm_gpu

    !> Computes hB = alpha*hA + beta*hB. Double precision
    !> @param[in] m     - The number of rows of the matrix A
    !> @param[in] n     - The number of columns of the matrix A
    !> @param[in] alpha - alpha scaling factor.
    !> @param[in] hA    - Host c-pointer to A matrix (entry).
    !> @param[in] lda   - The leading dimension of the array A.
    !> @param[in] beta  - beta scaling factor.
    !> @param[in,out] hB    - Host c-pointer to B matrix. On exit it have the result.
    !> @param[in] lda   - The leading dimension of the array A.
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine zgeadd_gpu(m, n, alpha, hA, lda, beta, hB, ldb, world)
        integer(i32), intent(in)            :: m
        integer(i32), intent(in)            :: n
        complex(r64), intent(in)            :: alpha
        type(c_ptr),  value                 :: hA
        integer(i32), intent(in)            :: lda
        complex(r64), intent(in)            :: beta
        type(c_ptr),  value                 :: hB
        integer(i32), intent(in)            :: ldb
        type(device_world_t), intent(inout) :: world

        complex(r64), pointer :: A(:,:), B(:,:)
        integer(i32) :: i, j

        call c_f_pointer(hA, A, [lda,n])
        call c_f_pointer(hB, B, [ldb,n])

#if defined(_INTEL_COMPILER)
        call mkl_zomatadd('c','n','n',m,n,A,lda,B,ldb,B,ldb)
#else
        !$omp parallel do collapse(2) default(shared) private(i,j)
        do i = 1, n
            do j = 1, m
                B(j,i) = beta * B(j,i) + alpha * A(j,i)
            end do
        end do
        !$omp end parallel do
#endif

        nullify(A,B)

    end subroutine zgeadd_gpu

    !> Copies x contents to y.
    !> @param[in]       n           - Number of elements in vectors x and y. n >= 0.
    !> @param[in]       hx          - Host pointer to x. The n element vector x of dimension (1 + (n-1)*incx).
    !> @param[in]       incx        - Stride between consecutive elements of hx.
    !> @param[in,out]   hy          - Host pointer to y. The n element vector y of dimension (1 + (n-1)*incy).
    !> @param[in]       incy        - Stride between consecutive elements of hy.
    !> @param[in,out]   world       - the device-host handler (CPU-backend)
    subroutine zcopy_gpu(n, hx, incx, hy, incy, world)
        integer(i32), intent(in)                        :: n
        type(c_ptr),  value                             :: hx
        integer(i32), intent(in)                        :: incx
        type(c_ptr),  value                             :: hy
        integer(i32), intent(in)                        :: incy
        type(device_world_t), intent(inout)             :: world

        complex(r64), pointer :: x(:), y(:)

        call c_f_pointer(hx, x, [1 + (n-1)*incx])
        call c_f_pointer(hy, y, [1 + (n-1)*incy])
        call zcopy(n, x, incx, y, incy)
        nullify(x,y)

    end subroutine zcopy_gpu

    !> Sets an array to zero
    !> @param[in,out] hA    - Host pointer to A array.
    !> @param[in,out] size  - the size of the array
    !> @param[in,out] world - the device-host handler. CPU backend.
    subroutine zsetzero_gpu(hA, size, world)
        type(c_ptr),  value                 :: hA
        integer(c_size_t), intent(in)       :: size
        type(device_world_t), intent(inout) :: world

        call memset(hA, 0_c_int, size*bytes_double_complex)

    end subroutine zsetzero_gpu

end module device_linalg_common_interface
