!
!   -- MAGMA (version 2.0) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date
!

module magma2

use iso_c_binding

use magma2_common
use magma2_sfortran
use magma2_dfortran
use magma2_cfortran
use magma2_zfortran

implicit none

!! =============================================================================
!! Parameter constants from magma_types.h
integer(c_int), parameter ::   &
    MagmaFalse         = 0,    &
    MagmaTrue          = 1,    &
    
    MagmaRowMajor      = 101,  &
    MagmaColMajor      = 102,  &
    
    MagmaNoTrans       = 111,  &
    MagmaTrans         = 112,  &
    MagmaConjTrans     = 113,  &
    
    MagmaUpper         = 121,  &
    MagmaLower         = 122,  &
    MagmaGeneral       = 123,  &
    MagmaFull          = 123,  &  !! deprecated, use MagmaGeneral
    
    MagmaNonUnit       = 131,  &
    MagmaUnit          = 132,  &
    
    MagmaLeft          = 141,  &
    MagmaRight         = 142,  &
    MagmaBothSides     = 143,  &

    MagmaHybrid        = 701,  &
    MagmaNative        = 702
!! todo all the rest


!! =============================================================================
!! Fortran interfaces to C functions
interface

    !! -------------------------------------------------------------------------
    !! initialize
    subroutine magma_init() &
    bind(C, name="magma_init")
        import
    end subroutine

    subroutine magma_finalize() &
    bind(C, name="magma_finalize")
        import
    end subroutine

    !! -------------------------------------------------------------------------
    !! version
    subroutine magma_version( major, minor, micro ) &
    bind(C, name="magma_version")
        import
        integer(c_int), target :: major, minor, micro
    end subroutine

    subroutine magma_print_environment() &
    bind(C, name="magma_print_environment")
        import
    end subroutine

    !! -------------------------------------------------------------------------
    !! timing
    real(c_double) function magma_wtime() &
    bind(C, name="magma_wtime")
        import
    end function

    real(c_double) function magma_sync_wtime( queue ) &
    bind(C, name="magma_wtime")
        import
        type(c_ptr), value :: queue
    end function

    !! -------------------------------------------------------------------------
    !! device support
    integer(c_int) function magma_num_gpus() &
    bind(C, name="magma_num_gpus")
        import
    end function

    integer(c_int) function magma_get_device_arch() &
    bind(C, name="magma_getdevice_arch")
        import
    end function

    subroutine magma_get_device( dev ) &
    bind(C, name="magma_getdevice")
        import
        integer(c_int), target :: dev
    end subroutine

    subroutine magma_set_device( dev ) &
    bind(C, name="magma_setdevice")
        import
        integer(c_int), value :: dev
    end subroutine

    integer(c_size_t) function magma_mem_size( queue ) &
    bind(C, name="magma_mem_size")
        import
        type(c_ptr), value :: queue
    end function

    !! -------------------------------------------------------------------------
    !! queue support
    subroutine magma_queue_create_internal( dev, queue_ptr, func, file, line ) &
    bind(C, name="magma_queue_create_internal")
        import
        integer(c_int), value :: dev
        type(c_ptr), target :: queue_ptr  !! queue_t*
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine

    subroutine magma_queue_destroy_internal( queue, func, file, line ) &
    bind(C, name="magma_queue_destroy_internal")
        import
        type(c_ptr), value :: queue  !! queue_t
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine

    subroutine magma_queue_sync_internal( queue, func, file, line ) &
    bind(C, name="magma_queue_sync_internal")
        import
        type(c_ptr), value :: queue  !! queue_t
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine

    integer(c_int) function magma_queue_get_device( queue ) &
    bind(C, name="magma_queue_get_device")
        import
        type(c_ptr), value :: queue  !! queue_t
    end function

    type(c_ptr) function magma_queue_get_cuda_stream( queue ) &
    bind(C, name="magma_queue_get_cuda_stream")
        import
        type(c_ptr), value :: queue  !! queue_t
    end function

    type(c_ptr) function magma_queue_get_hip_stream( queue ) &
    bind(C, name="magma_queue_get_hip_stream")
        import
        type(c_ptr), value :: queue  !! queue_t
    end function

    !! -------------------------------------------------------------------------
    !! offsets pointers -- 1D vectors with inc
    !! see offset.c
    type(c_ptr) function magma_soffset_1d( ptr, inc, i ) &
    bind(C, name="magma_soffset_1d")
        import
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_doffset_1d( ptr, inc, i ) &
    bind(C, name="magma_doffset_1d")
        import
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_coffset_1d( ptr, inc, i ) &
    bind(C, name="magma_coffset_1d")
        import
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_zoffset_1d( ptr, inc, i ) &
    bind(C, name="magma_zoffset_1d")
        import
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_ioffset_1d( ptr, inc, i ) &
    bind(C, name="magma_ioffset_1d")
        import
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    !! -------------------------------------------------------------------------
    !! offsets pointers -- 2D matrices with lda
    !! see offset.c
    type(c_ptr) function magma_soffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_soffset_2d")
        import
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_doffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_doffset_2d")
        import
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_coffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_coffset_2d")
        import
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_zoffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_zoffset_2d")
        import
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_ioffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_ioffset_2d")
        import
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

end interface

!! =============================================================================
!! Fortran routines & functions
contains

    !! -------------------------------------------------------------------------
    !! queue support
    subroutine magma_queue_create( dev, queue_ptr )
        integer(c_int), value :: dev
        type(c_ptr), target :: queue_ptr  !! queue_t*
        
        call magma_queue_create_internal( &
                dev, queue_ptr, &
                "magma_queue_create" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_queue_destroy( queue )
        type(c_ptr), value :: queue  !! queue_t
        
        call magma_queue_destroy_internal( &
                queue, &
                "magma_queue_destroy" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_queue_sync( queue )
        type(c_ptr), value :: queue  !! queue_t
        
        call magma_queue_sync_internal( &
                queue, &
                "magma_queue_sync" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

end module
