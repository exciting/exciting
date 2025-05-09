!! @generated from magma2_zfortran.F90, fortran z -> c, Wed Nov  1 17:17:00 2023

module magma2_cfortran

use iso_c_binding
use magma2_common

implicit none

!! =============================================================================
!! Fortran interfaces to C functions
interface

    !! -------------------------------------------------------------------------
    !! CPU interfaces (matrix in CPU memory)
    subroutine magma_cgetrf( m, n, A, lda, ipiv, info ) &
    bind(C, name="magma_cgetrf")
        import
        integer(c_int),            value  :: m, n, lda
        complex(c_float_complex), target :: A(lda,*)
        integer(c_int),            target :: ipiv(*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    subroutine magma_cpotrf( uplo, n, A, lda, info ) &
    bind(C, name="magma_cpotrf")
        import
        integer(c_int),            value  :: uplo
        integer(c_int),            value  :: n, lda
        complex(c_float_complex), target :: A(lda,*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    !! -------------------------------------------------------------------------
    !! GPU interfaces (matrix in GPU memory)
    subroutine magma_cgetrf_gpu( m, n, dA, lda, ipiv, info ) &
    bind(C, name="magma_cgetrf_gpu")
        import
        integer(c_int), value  :: m, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: ipiv(*)
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_cpotrf_gpu( uplo, n, dA, lda, info ) &
    bind(C, name="magma_cpotrf_gpu")
        import
        integer(c_int), value  :: uplo, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_cgetri_gpu( n, dA, lda, ipiv, dwork, lwork, info ) &
        bind(C, name="magma_cgetri_gpu")
            import
            integer(c_int), value  :: n
            type(c_ptr),    value  :: dA
            integer(c_int), value  :: lda
            integer(c_int), target :: ipiv(*)
            type(c_ptr),    value  :: dwork
            integer(c_int), value  :: lwork
            integer(c_int), target :: info  !! int*
    end subroutine

    function magma_get_cgetri_nb( n ) result(size) &
        bind(C, name="magma_get_cgetri_nb")
            import
            integer(c_int) :: n
            integer(c_int) :: size
    end function 

    !! -------------------------------------------------------------------------
    !! batched GPU interfaces (all arrays in GPU memory)
    subroutine magma_cgetrf_batched( &
        m, n, dA_array, lda, ipiv_array, info_array, batchcount, queue ) &
    bind(C, name="magma_cgetrf_batched")
        import
        integer(c_int), value  :: m, n, lda, batchcount
        type(c_ptr),    value  :: dA_array    !! double_complex**
        type(c_ptr),    value  :: ipiv_array  !! int**
        type(c_ptr),    value  :: info_array  !! int*
        type(c_ptr),    value  :: queue
    end subroutine

    subroutine magma_cgemm_batched_strided(transA, transB, m, n, k, alpha, dA, lda, strideA, &
                                   dB, ldb, strideB, beta, dC, ldc, strideC, batchCount, queue) &
    bind(C, name="magmablas_cgemm_batched_strided")
        import
        integer(c_int),             value :: transA, transB, m, n, k, lda, ldb, ldc, batchCount
        integer(c_int),             value :: strideA, strideB, strideC
        complex(c_float_complex),   value :: alpha, beta
        type(c_ptr),                value :: dA, dB, dC
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    !! -------------------------------------------------------------------------
    !! BLAS (matrices in GPU memory)
    subroutine magma_caxpy( &
        n, &
        alpha, dx, incx, &
               dy, incy, &
        queue ) &
    bind(C, name="magma_caxpy")
        import
        integer(c_int),             value :: n, incx, incy
        complex(c_float_complex),   value :: alpha
        type(c_ptr),                value :: dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_cgeru(m, n, alpha, dx, incx, dy, incy, da, ldda, queue ) &
    bind(C, name="magma_cgeru")
        import
        integer(c_int),             value :: m, n, incx, incy, ldda
        complex(c_float_complex),   value :: alpha
        type(c_ptr),                value :: dx, dy, da
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_cgerc(m, n, alpha, dx, incx, dy, incy, da, ldda, queue ) &
        bind(C, name="magma_cgerc")
            import
            integer(c_int),             value :: m, n, incx, incy, ldda
            complex(c_float_complex),   value :: alpha
            type(c_ptr),                value :: dx, dy, da
            type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    complex(c_float_complex) function magma_cdotc( &
        n, dx, incx, dy, incy, queue) bind(C, name="magma_cdotc")
        import
        integer(c_int),             value :: n, incx, incy
        type(c_ptr),                value :: dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end function magma_cdotc

    complex(c_float_complex) function magma_cdotu( &
        n, dx, incx, dy, incy, queue) bind(C, name="magma_cdotu")
        import
        integer(c_int),             value :: n, incx, incy
        type(c_ptr),                value :: dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end function magma_cdotu

    subroutine magma_cgemv( &
        transA, m, n, &
        alpha, dA, lda, &
               dx, incx, &
        beta,  dy, incy, &
        queue ) &
    bind(C, name="magma_cgemv")
        import
        integer(c_int),             value :: transA, m, n, lda, incx, incy
        complex(c_float_complex),  value :: alpha, beta
        type(c_ptr),                value :: dA, dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_cgemm( &
        transA, transB, m, n, k, &
        alpha, dA, lda, &
               dB, ldb, &
        beta,  dC, ldc, &
        queue ) &
    bind(C, name="magma_cgemm")
        import
        integer(c_int),             value :: transA, transB, m, n, k, lda, ldb, ldc
        complex(c_float_complex),  value :: alpha, beta
        type(c_ptr),                value :: dA, dB, dC
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_chemm( &
        side, uplo, m, n,   &
        alpha, dA, lda, &
               dB, ldb, &
        beta,  dC, ldc, &
        queue ) &
    bind(C, name="magma_chemm")
        import
        integer(c_int),             value :: side, uplo, m, n, lda, ldb, ldc
        complex(c_float_complex),   value :: alpha, beta
        type(c_ptr),                value :: dA, dB, dC
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magmablas_cgeadd2(m, n, alpha, dA, lda, beta, dB, ldb, queue) &
    bind(C, name="magmablas_cgeadd2")
        import
        integer(c_int),             value :: m, n, lda, ldb
        complex(c_float_complex),   value :: alpha, beta
        type(c_ptr),                value :: dA, dB
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_chemv(uplo, n, alpha, dA, lda, dx, incx, beta, dy, incy, queue ) &
    bind(C, name="magma_chemv")
        import
        integer(c_int),             value :: uplo, n, lda, incx, incy
        complex(c_float_complex),   value :: alpha, beta
        type(c_ptr),                value :: dA, dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

end interface

!! =============================================================================
!! Fortran routines & functions
contains

    !! -------------------------------------------------------------------------
    !! malloc wrappers
    integer(c_int) function magma_cmalloc( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_cmalloc = magma_malloc( ptr, n*sizeof_complex )
    end function

    integer(c_int) function magma_cmalloc_cpu( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_cmalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_complex )
    end function

    integer(c_int) function magma_cmalloc_pinned( ptr, n )

        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_cmalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_complex )
    end function

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    subroutine magma_csetmatrix( &
        m, n, hA_src, lda, dB_dst, ldb, queue )

        integer(c_int),            value  :: m, n, lda, ldb
        complex(c_float_complex), target :: hA_src(lda,*)
        type(c_ptr),               value  :: dB_dst
        type(c_ptr),               value  :: queue
        
        call magma_setmatrix_internal( &
                m, n, int(sizeof_complex), c_loc(hA_src), lda, dB_dst, ldb, queue, &
                "magma_csetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_cgetmatrix( &
        m, n, dA_src, lda, hB_dst, ldb, queue )

        integer(c_int),            value  :: m, n, lda, ldb
        type(c_ptr),               value  :: dA_src
        complex(c_float_complex), target :: hB_dst(ldb,*)
        type(c_ptr),               value  :: queue
        
        call magma_getmatrix_internal( &
                m, n, int(sizeof_complex), dA_src, lda, c_loc(hB_dst), ldb, queue, &
                "magma_cgetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    !! -----------------------------------------------------------------------------
    !! Copy
    subroutine magma_ccopyvector_async(n, dx_src, incx, dy_dst, incy, queue)
        integer(c_int), value :: n, incx, incy
        type(c_ptr), value :: dx_src, dy_dst, queue

        call magma_copyvector_async_internal(n, int(sizeof_complex), dx_src, incx, &
                                             dy_dst, incy, queue, &
                                             "magma_ccopyvector_async" // c_null_char, &
                                              __FILE__ // c_null_char, &
                                              __LINE__ )

    end subroutine magma_ccopyvector_async

end module
