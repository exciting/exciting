!> Module providing a unified interface for FFT depending
!> on the compilation choices.
module m_zfftifc

            use precision, only: i32, dp
#if FFTW
            use, intrinsic :: iso_c_binding
            include 'fftw3.f03'
#else
            implicit none
#endif

            private
            public :: zfftifc

contains


      subroutine zfftifc (nd, n, sgn, z)

      !> number of dimensions
      integer(i32), intent(in) :: nd
      !> FFT direction, -1: forward; 1: backward
      integer(i32), intent(in) :: sgn
      !> grid sizes
      integer(i32), intent(in) :: n (nd)
      !> array to transform
      complex (dp), intent(inout) :: z (*)

      !-------------------------------------!
      !     interface to FFTW version 3     !
      !-------------------------------------!

#ifdef FFTW

      integer(i32) :: i
      type(c_ptr)  :: plan
      real(dp)     :: scaling
      ! FFTW3 init globals
      integer(c_int) :: cerror

      ! Note that for vanilla FFTW3 plan creation and destruction are not
      ! thread safe and thus need to be protected by a critical statement
      ! See https://www.fftw.org/fftw3_doc/Thread-safety.html
      !$omp critical
      call dfftw_plan_dft(plan, nd, n, z, z, sgn, FFTW_ESTIMATE)
      !$omp end critical

      call dfftw_execute(plan)

      !$omp critical
      call dfftw_destroy_plan(plan)
      !$omp end critical

      ! FFTW does not normalize the result of the inverse transform;
      ! manual scaling is required
      if (sgn == -1) then
            scaling = 1.0_dp / real(product(n(1:nd)), kind=dp)
            do i = 1, product(n(1:nd))
                  z(i) = scaling * z(i)
            end do
      end if

#else
      !----------------------------------------!
      !     interface to modified FFTPACK5     !
      !----------------------------------------!
            Call cfftnd (nd, n, sgn, z)
#endif

      end subroutine zfftifc

end module m_zfftifc
