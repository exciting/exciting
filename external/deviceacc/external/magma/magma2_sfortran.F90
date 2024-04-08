!! @generated from magma2_zfortran.F90, fortran z -> s, Wed Nov  1 17:17:00 2023

module magma2_sfortran

use iso_c_binding
use magma2_common
implicit none

!! =============================================================================
!! Fortran interfaces to C functions
interface

    !! -------------------------------------------------------------------------
    !! CPU interfaces (matrix in CPU memory)
    subroutine magma_sgetrf( m, n, A, lda, ipiv, info ) &
    bind(C, name="magma_sgetrf")
        import
        integer(c_int),            value  :: m, n, lda
        real(c_float), target :: A(lda,*)
        integer(c_int),            target :: ipiv(*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    subroutine magma_spotrf( uplo, n, A, lda, info ) &
    bind(C, name="magma_spotrf")
        import
        integer(c_int),            value  :: uplo
        integer(c_int),            value  :: n, lda
        real(c_float), target :: A(lda,*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    !! -------------------------------------------------------------------------
    !! GPU interfaces (matrix in GPU memory)
    subroutine magma_sgetrf_gpu( m, n, dA, lda, ipiv, info ) &
    bind(C, name="magma_sgetrf_gpu")
        import
        integer(c_int), value  :: m, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: ipiv(*)
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_spotrf_gpu( uplo, n, dA, lda, info ) &
    bind(C, name="magma_spotrf_gpu")
        import
        integer(c_int), value  :: uplo, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: info  !! int*
    end subroutine

    !! -------------------------------------------------------------------------
    !! batched GPU interfaces (all arrays in GPU memory)
    subroutine magma_sgetrf_batched( &
        m, n, dA_array, lda, ipiv_array, info_array, batchcount, queue ) &
    bind(C, name="magma_sgetrf_batched")
        import
        integer(c_int), value  :: m, n, lda, batchcount
        type(c_ptr),    value  :: dA_array    !! double_real**
        type(c_ptr),    value  :: ipiv_array  !! int**
        type(c_ptr),    value  :: info_array  !! int*
        type(c_ptr),    value  :: queue
    end subroutine

    !! -------------------------------------------------------------------------
    !! BLAS (matrices in GPU memory)
    subroutine magma_saxpy( &
        n, &
        alpha, dx, incx, &
               dy, incy, &
        queue ) &
    bind(C, name="magma_saxpy")
        import
        integer(c_int),             value :: n, incx, incy
        real(c_float),  value :: alpha
        type(c_ptr),                value :: dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_sgemv( &
        transA, m, n, &
        alpha, dA, lda, &
               dx, incx, &
        beta,  dy, incy, &
        queue ) &
    bind(C, name="magma_sgemv")
        import
        integer(c_int),             value :: transA, m, n, lda, incx, incy
        real(c_float),  value :: alpha, beta
        type(c_ptr),                value :: dA, dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_sgemm( &
        transA, transB, m, n, k, &
        alpha, dA, lda, &
               dB, ldb, &
        beta,  dC, ldc, &
        queue ) &
    bind(C, name="magma_sgemm")
        import
        integer(c_int),             value :: transA, transB, m, n, k, lda, ldb, ldc
        real(c_float),  value :: alpha, beta
        type(c_ptr),                value :: dA, dB, dC
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

end interface

!! =============================================================================
!! Fortran routines & functions
contains

    !! -------------------------------------------------------------------------
    !! malloc wrappers
    integer(c_int) function magma_smalloc( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_smalloc = magma_malloc( ptr, n*sizeof_real )
    end function

    integer(c_int) function magma_smalloc_cpu( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_smalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_real )
    end function

    integer(c_int) function magma_smalloc_pinned( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_smalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_real )
    end function

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    subroutine magma_ssetmatrix( &
        m, n, hA_src, lda, dB_dst, ldb, queue )

        integer(c_int),            value  :: m, n, lda, ldb
        real(c_float), target :: hA_src(lda,*)
        type(c_ptr),               value  :: dB_dst
        type(c_ptr),               value  :: queue
        
        call magma_setmatrix_internal( &
                m, n, int(sizeof_real), c_loc(hA_src), lda, dB_dst, ldb, queue, &
                "magma_ssetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_sgetmatrix( &
        m, n, dA_src, lda, hB_dst, ldb, queue )

        integer(c_int),            value  :: m, n, lda, ldb
        type(c_ptr),               value  :: dA_src
        real(c_float), target :: hB_dst(ldb,*)
        type(c_ptr),               value  :: queue
        
        call magma_getmatrix_internal( &
                m, n, int(sizeof_real), dA_src, lda, c_loc(hB_dst), ldb, queue, &
                "magma_sgetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

end module
