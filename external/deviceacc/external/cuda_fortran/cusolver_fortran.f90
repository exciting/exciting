module cusolver_fortran

    use iso_c_binding
    use cublas_fortran

    implicit none

    enum, bind(c)
        enumerator :: CUSOLVER_STATUS_SUCCESS                   = 0
        enumerator :: CUSOLVER_STATUS_NOT_INITIALIZED           = 1
        enumerator :: CUSOLVER_STATUS_ALLOC_FAILED              = 2
        enumerator :: CUSOLVER_STATUS_INVALID_VALUE             = 3
        enumerator :: CUSOLVER_STATUS_ARCH_MISMATCH             = 4
        enumerator :: CUSOLVER_STATUS_MAPPING_ERROR             = 5
        enumerator :: CUSOLVER_STATUS_EXECUTION_FAILED          = 6
        enumerator :: CUSOLVER_STATUS_INTERNAL_ERROR            = 7
        enumerator :: CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED = 8
        enumerator :: CUSOLVER_STATUS_NOT_SUPPORTED             = 9
        enumerator :: CUSOLVER_STATUS_ZERO_PIVOT                = 10
        enumerator :: CUSOLVER_STATUS_INVALID_LICENSE           = 11
        enumerator :: CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED= 12
        enumerator :: CUSOLVER_STATUS_IRS_PARAMS_INVALID        = 13
        enumerator :: CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC   = 14
        enumerator :: CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE = 15
        enumerator :: CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER= 16
        enumerator :: CUSOLVER_STATUS_IRS_INTERNAL_ERROR        = 20
        enumerator :: CUSOLVER_STATUS_IRS_NOT_SUPPORTED         = 21
        enumerator :: CUSOLVER_STATUS_IRS_OUT_OF_RANGE          = 22
        enumerator :: CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES=23
        enumerator :: CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED = 25
        enumerator :: CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED   = 26
        enumerator :: CUSOLVER_STATUS_IRS_MATRIX_SINGULAR       = 30
        enumerator :: CUSOLVER_STATUS_INVALID_WORKSPACE         = 31
    end enum

    enum, bind(c)
        enumerator :: CUSOLVER_EIG_TYPE_1 = 1
        enumerator :: CUSOLVER_EIG_TYPE_2 = 2
        enumerator :: CUSOLVER_EIG_TYPE_3 = 3
    end enum

    enum, bind(c)
        enumerator :: CUSOLVER_EIG_MODE_NOVECTOR = 0
        enumerator :: CUSOLVER_EIG_MODE_VECTOR   = 1
    end enum

    enum, bind(c)
        enumerator :: CUSOLVER_EIG_RANGE_ALL = 1001
        enumerator :: CUSOLVER_EIG_RANGE_I   = 1002
        enumerator :: CUSOLVER_EIG_RANGE_V   = 1003
    end enum

    enum, bind(c)
        enumerator :: CUSOLVER_INF_NORM = 104
        enumerator :: CUSOLVER_MAX_NORM = 105
        enumerator :: CUSOLVER_ONE_NORM = 106
        enumerator :: CUSOLVER_FRO_NORM = 107
    end enum

    enum, bind(c)
        enumerator :: CUSOLVER_IRS_REFINE_NOT_SET          = 1100
        enumerator :: CUSOLVER_IRS_REFINE_NONE             = 1101
        enumerator :: CUSOLVER_IRS_REFINE_CLASSICAL        = 1102
        enumerator :: CUSOLVER_IRS_REFINE_CLASSICAL_GMRES  = 1103
        enumerator :: CUSOLVER_IRS_REFINE_GMRES            = 1104
        enumerator :: CUSOLVER_IRS_REFINE_GMRES_GMRES      = 1105
        enumerator :: CUSOLVER_IRS_REFINE_GMRES_NOPCOND    = 1106

        enumerator :: CUSOLVER_PREC_DD           = 1150
        enumerator :: CUSOLVER_PREC_SS           = 1151
        enumerator :: CUSOLVER_PREC_SHT          = 1152
    end enum

    enum, bind(c)
        enumerator :: CUSOLVER_R_8I  = 1201
        enumerator :: CUSOLVER_R_8U  = 1202
        enumerator :: CUSOLVER_R_64F = 1203
        enumerator :: CUSOLVER_R_32F = 1204
        enumerator :: CUSOLVER_R_16F = 1205
        enumerator :: CUSOLVER_R_16BF  = 1206
        enumerator :: CUSOLVER_R_TF32  = 1207
        enumerator :: CUSOLVER_R_AP  = 1208
        enumerator :: CUSOLVER_C_8I  = 1211
        enumerator :: CUSOLVER_C_8U  = 1212
        enumerator :: CUSOLVER_C_64F = 1213
        enumerator :: CUSOLVER_C_32F = 1214
        enumerator :: CUSOLVER_C_16F = 1215
        enumerator :: CUSOLVER_C_16BF  = 1216
        enumerator :: CUSOLVER_C_TF32  = 1217
        enumerator :: CUSOLVER_C_AP  = 1218
    end enum

    enum, bind(c)
        enumerator :: CUSOLVER_ALG_0 = 0
        enumerator :: CUSOLVER_ALG_1 = 1
    end enum

    interface cusolverDnCreate
        function cusolverDnCreate_(handle) bind(c, name="cusolverDnCreate")
            import :: c_ptr, c_int
            integer(c_int) :: cusolverDnCreate_
            type(c_ptr), value :: handle
        end function cusolverDnCreate_
    end interface cusolverDnCreate

    interface cusolverDnDestroy
        function cusolverDnDestroy_(handle) bind(c, name="cusolverDnDestroy")
            import :: c_ptr, c_int
            integer(c_int) :: cusolverDnDestroy_
            type(c_ptr), value :: handle
        end function cusolverDnDestroy_
    end interface cusolverDnDestroy

    interface cusolverDnChegvdx_bufferSize
        integer(c_int) function cusolverDnChegvdx_bufferSize_(handle, itype, jobz, range, uplo, &
                                            n, A, lda, B, ldb, vl, vu, il, iu, &
                                            h_meig, W, lwork) bind(c, name="cusolverDnChegvdx_bufferSize")
            import :: c_ptr, c_int, c_float
            type(c_ptr),    value       :: handle       ! cusolverDnHandle_t
            integer(c_int), intent(in)  :: itype        ! cusolverEigType_t
            integer(c_int), intent(in)  :: jobz         ! cusolverEigMode_t
            integer(c_int), intent(in)  :: range        ! cusolverEigRange_t
            integer(c_int), intent(in)  :: uplo         ! cublasFillMode_t
            integer(c_int), intent(in)  :: n
            type(c_ptr),    value       :: A  ! Device
            integer(c_int), intent(in)  :: lda
            type(c_ptr),    value       :: B  ! Device
            integer(c_int), intent(in)  :: ldb
            real(c_float), intent(in)   :: vl
            real(c_float), intent(in)   :: vu
            integer(c_int), intent(in)  :: il
            integer(c_int), intent(in)  :: iu
            integer(c_int), intent(out) :: h_meig           ! output
            type(c_ptr),    value       :: W                ! Device
            integer(c_int), intent(out) :: lwork            ! output
        end function cusolverDnChegvdx_bufferSize_
    end interface cusolverDnChegvdx_bufferSize

    interface cusolverDnChegvdx
        integer(c_int) function cusolverDnChegvdx_(handle, itype, jobz, range, uplo, &
                                n, A, lda, B, ldb, vl, vu, il, iu, &
                                h_meig, W, work, lwork, devInfo) bind(c, name="cusolverDnChegvdx")
            import :: c_ptr, c_int, c_float
            type(c_ptr), value           :: handle        ! cusolverDnHandle_t
            integer(c_int), intent(in)   :: itype         ! cusolverEigType_t
            integer(c_int), intent(in)   :: jobz          ! cusolverEigMode_t
            integer(c_int), intent(in)   :: range         ! cusolverEigRange_t
            integer(c_int), intent(in)   :: uplo          ! cublasFillMode_t
            integer(c_int), intent(in)   :: n
            type(c_ptr), value           :: A  ! Device
            integer(c_int), intent(in)   :: lda
            type(c_ptr), value           :: B  ! Device
            integer(c_int), intent(in)   :: ldb
            real(c_float), intent(in)    :: vl
            real(c_float), intent(in)    :: vu
            integer(c_int), intent(in)   :: il
            integer(c_int), intent(in)   :: iu
            integer(c_int), intent(out)  :: h_meig     ! output
            type(c_ptr), value           :: W          ! Device
            type(c_ptr), value           :: work       ! cuDoubleComplex*
            integer(c_int), intent(in)   :: lwork
            integer(c_int), intent(out)  :: devInfo    ! output (Device)
        end function cusolverDnChegvdx_
    end interface cusolverDnChegvdx

    interface cusolverDnZhegvdx_bufferSize
        integer(c_int) function cusolverDnZhegvdx_bufferSize_(handle, itype, jobz, range, uplo, &
                                            n, A, lda, B, ldb, vl, vu, il, iu, &
                                            h_meig, W, lwork) bind(c, name="cusolverDnZhegvdx_bufferSize")
            import :: c_ptr, c_int, c_double
            type(c_ptr),    value       :: handle       ! cusolverDnHandle_t
            integer(c_int), intent(in)  :: itype        ! cusolverEigType_t
            integer(c_int), intent(in)  :: jobz         ! cusolverEigMode_t
            integer(c_int), intent(in)  :: range        ! cusolverEigRange_t
            integer(c_int), intent(in)  :: uplo         ! cublasFillMode_t
            integer(c_int), intent(in)  :: n
            type(c_ptr),    value       :: A  ! Device
            integer(c_int), intent(in)  :: lda
            type(c_ptr),    value       :: B  ! Device
            integer(c_int), intent(in)  :: ldb
            real(c_double), intent(in)  :: vl
            real(c_double), intent(in)  :: vu
            integer(c_int), intent(in)  :: il
            integer(c_int), intent(in)  :: iu
            integer(c_int), intent(out) :: h_meig           ! output
            type(c_ptr),    value       :: W                ! Device
            integer(c_int), intent(out) :: lwork            ! output
        end function cusolverDnZhegvdx_bufferSize_
    end interface cusolverDnZhegvdx_bufferSize

    interface cusolverDnZhegvdx
        integer(c_int) function cusolverDnZhegvdx_(handle, itype, jobz, range, uplo, &
                                n, A, lda, B, ldb, vl, vu, il, iu, &
                                h_meig, W, work, lwork, devInfo) bind(c, name="cusolverDnZhegvdx")
            import :: c_ptr, c_int, c_double
            type(c_ptr), value           :: handle        ! cusolverDnHandle_t
            integer(c_int), intent(in)   :: itype         ! cusolverEigType_t
            integer(c_int), intent(in)   :: jobz          ! cusolverEigMode_t
            integer(c_int), intent(in)   :: range         ! cusolverEigRange_t
            integer(c_int), intent(in)   :: uplo          ! cublasFillMode_t
            integer(c_int), intent(in)   :: n
            type(c_ptr), value           :: A  ! Device
            integer(c_int), intent(in)   :: lda
            type(c_ptr), value           :: B  ! Device
            integer(c_int), intent(in)   :: ldb
            real(c_double), intent(in)   :: vl
            real(c_double), intent(in)   :: vu
            integer(c_int), intent(in)   :: il
            integer(c_int), intent(in)   :: iu
            integer(c_int), intent(out)  :: h_meig     ! output
            type(c_ptr), value           :: W          ! Device
            type(c_ptr), value           :: work       ! cuDoubleComplex*
            integer(c_int), intent(in)   :: lwork
            integer(c_int), intent(out)  :: devInfo    ! output (Device)
        end function cusolverDnZhegvdx_
    end interface cusolverDnZhegvdx

end module cusolver_fortran
