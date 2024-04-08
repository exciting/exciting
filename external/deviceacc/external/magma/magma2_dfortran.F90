!! @generated from magma2_zfortran.F90, fortran z -> d, Wed Nov  1 17:17:00 2023

module magma2_dfortran

use iso_c_binding
use magma2_common

implicit none

!! =============================================================================
!! Fortran interfaces to C functions
interface

    !! -------------------------------------------------------------------------
    !! CPU interfaces (matrix in CPU memory)
    subroutine magma_dgetrf( m, n, A, lda, ipiv, info ) &
    bind(C, name="magma_dgetrf")
        import
        integer(c_int),            value  :: m, n, lda
        real(c_double), target :: A(lda,*)
        integer(c_int),            target :: ipiv(*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    subroutine magma_dpotrf( uplo, n, A, lda, info ) &
    bind(C, name="magma_dpotrf")
        import
        integer(c_int),            value  :: uplo
        integer(c_int),            value  :: n, lda
        real(c_double), target :: A(lda,*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    !! -------------------------------------------------------------------------
    !! GPU interfaces (matrix in GPU memory)
    subroutine magma_dgetrf_gpu( m, n, dA, lda, ipiv, info ) &
    bind(C, name="magma_dgetrf_gpu")
        import
        integer(c_int), value  :: m, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: ipiv(*)
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_dpotrf_gpu( uplo, n, dA, lda, info ) &
    bind(C, name="magma_dpotrf_gpu")
        import
        integer(c_int), value  :: uplo, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: info  !! int*
    end subroutine

    !! -------------------------------------------------------------------------
    !! batched GPU interfaces (all arrays in GPU memory)
    subroutine magma_dgetrf_batched( &
        m, n, dA_array, lda, ipiv_array, info_array, batchcount, queue ) &
    bind(C, name="magma_dgetrf_batched")
        import
        integer(c_int), value  :: m, n, lda, batchcount
        type(c_ptr),    value  :: dA_array    !! double_real**
        type(c_ptr),    value  :: ipiv_array  !! int**
        type(c_ptr),    value  :: info_array  !! int*
        type(c_ptr),    value  :: queue
    end subroutine

    !! -------------------------------------------------------------------------
    !! BLAS (matrices in GPU memory)
    subroutine magma_daxpy( &
        n, &
        alpha, dx, incx, &
               dy, incy, &
        queue ) &
    bind(C, name="magma_daxpy")
        import
        integer(c_int),             value :: n, incx, incy
        real(c_double),  value :: alpha
        type(c_ptr),                value :: dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_dgemv( &
        transA, m, n, &
        alpha, dA, lda, &
               dx, incx, &
        beta,  dy, incy, &
        queue ) &
    bind(C, name="magma_dgemv")
        import
        integer(c_int),             value :: transA, m, n, lda, incx, incy
        real(c_double),  value :: alpha, beta
        type(c_ptr),                value :: dA, dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_dgemm( &
        transA, transB, m, n, k, &
        alpha, dA, lda, &
               dB, ldb, &
        beta,  dC, ldc, &
        queue ) &
    bind(C, name="magma_dgemm")
        import
        integer(c_int),             value :: transA, transB, m, n, k, lda, ldb, ldc
        real(c_double),  value :: alpha, beta
        type(c_ptr),                value :: dA, dB, dC
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

end interface

!! =============================================================================
!! Fortran routines & functions
contains

    !! -------------------------------------------------------------------------
    !! malloc wrappers
    integer(c_int) function magma_dmalloc( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_dmalloc = magma_malloc( ptr, n*sizeof_double )
    end function

    integer(c_int) function magma_dmalloc_cpu( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_dmalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_double )
    end function

    integer(c_int) function magma_dmalloc_pinned( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_dmalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_double )
    end function

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    subroutine magma_dsetmatrix( &
        m, n, hA_src, lda, dB_dst, ldb, queue )

        integer(c_int),            value  :: m, n, lda, ldb
        real(c_double), target :: hA_src(lda,*)
        type(c_ptr),               value  :: dB_dst
        type(c_ptr),               value  :: queue
        
        call magma_setmatrix_internal( &
                m, n, int(sizeof_double), c_loc(hA_src), lda, dB_dst, ldb, queue, &
                "magma_dsetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_dgetmatrix( &
        m, n, dA_src, lda, hB_dst, ldb, queue )

        integer(c_int),            value  :: m, n, lda, ldb
        type(c_ptr),               value  :: dA_src
        real(c_double), target :: hB_dst(ldb,*)
        type(c_ptr),               value  :: queue
        
        call magma_getmatrix_internal( &
                m, n, int(sizeof_double), dA_src, lda, c_loc(hB_dst), ldb, queue, &
                "magma_dgetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

end module
