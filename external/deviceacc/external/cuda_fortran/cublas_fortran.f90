module cublas_fortran

    use iso_c_binding

    implicit none

    enum, bind(c)
        enumerator :: CUBLAS_STOREV_COLUMNWISE = 0
        enumerator :: CUBLAS_STOREV_ROWWISE    = 1
    end enum

    enum, bind(c)
        enumerator :: CUBLAS_DIRECT_FORWARD    = 0
        enumerator :: CUBLAS_DIRECT_BACKWARD   = 1
    end enum

    enum, bind(c)
        enumerator :: CUBLAS_FILL_MODE_LOWER=0
        enumerator :: CUBLAS_FILL_MODE_UPPER=1
    end enum

end module cublas_fortran
