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
!> This file contains a type to perform FFT using device accelerated
!> routines. The user does not need to take care of the device vendor
!> Supported vendors include NVIDIA, and AMD. This is the host backend.
module m_fft_device
    
    use iso_c_binding,    only: c_int, c_ptr, c_null_ptr, c_double_complex, c_float_complex, &
                                c_int32_t, c_intptr_t, c_double, c_float, c_funptr, c_size_t, &
                                c_f_pointer, c_char
                        
    use iso_fortran_env,  only: i32=>int32, r32=>real32, r64=>real64
    use m_device_world_t, only: device_world_t

    implicit none

    include "fftw3.f03"

    private 
    public :: fft_device_t

    !> Type for setting up FFT plans.
    type fft_device_t
        !> Dimensions of the FFT.
        integer(c_int), allocatable, private :: dims(:)
        !> Rank.
        integer, private :: rank = 0
        !> Sign of the FFT. The definition depends on the backend library
        integer, private :: fft_sign = 0
        !> Actual FFT plan.
        type(c_ptr), private :: plan = c_null_ptr
        !> Double precision
        logical, private :: is_double
        !> Device world
        type(device_world_t), pointer :: world
    contains
        procedure :: initialize
        procedure :: delete
        procedure :: execute
    end type fft_device_t

contains 

    !> Initializes the FFT type using device accelerated procedures
    !> @param[out] this    - the fft_device_t to initialize
    !> @param[in]  dims    - the shape of the problem to solve
    !> @param[in]  fft_sign    - the direction of the Fourier transform (forward: -1 // backward: 1) 
    !> @param[in,out]  df      - the data or that of similar shape over which the FFT is computed (host pointer)
    !> @param[in]  is_double  - are we dealing with double precision, be aware as we use c_ptr the type information is lost
    !> @param[in,out] world - device handler 
    subroutine initialize(this, dims, fft_sign, df, is_double, world)
        
        class(fft_device_t), intent(out)            :: this
        integer, intent(in)                         :: dims(:)
        integer, intent(in)                         :: fft_sign
        type(c_ptr), intent(inout)                  :: df
        logical, intent(in)                         :: is_double 
        type(device_world_t), target, intent(inout) :: world 

        complex(c_double_complex), pointer :: fortran_df_double(:)
        complex(c_float_complex),  pointer :: fortran_df_single(:)

        this%rank   = int(size(dims), kind = c_int)
        this%dims   = int(dims(this%rank:1:-1), kind = c_int) ! Used libraries expect C-order
        this%fft_sign   = fft_sign
        this%is_double = is_double
        this%world  => world

        ! Some sanity check
        if (this%fft_sign /= -1 .and. this%fft_sign /= 1) then
            error stop "fft_type%initialize: The only accepted fft_signs are -1 (forward) and 1 (backward)."
        end if

        if (this%is_double) then
            call c_f_pointer(df, fortran_df_double, [product(this%dims)])
            this%plan = fftw_plan_dft( &
                        rank=this%rank, &
                        n=this%dims, &
                        in=fortran_df_double, &
                        out=fortran_df_double, &
                        sign=this%fft_sign, &
                        flags=FFTW_ESTIMATE)
            nullify(fortran_df_double)
        else
            call c_f_pointer(df, fortran_df_single, [product(this%dims)])
            this%plan = fftwf_plan_dft( &
                    rank=this%rank, &
                    n=this%dims, &
                    in=fortran_df_single, &
                    out=fortran_df_single, &
                    sign=this%fft_sign, &
                    flags=FFTW_ESTIMATE)
            nullify(fortran_df_single)
        end if

    end subroutine initialize

    !> Executes the plan
    !> @param[in,out]      this - the fft_device_t for which to execute the plan
    !> @param[in,out]  df   - the data over which to perform the plan (host ptr). In exit contains the result.
    !> @param[in]      rescale_forward - rescale the FFT in case of forward FFT
    !> @param[in]      synchronize - force syncronization after the execution of the plan
    subroutine execute(this, df, rescale_forward, synchronize)

        class(fft_device_t), intent(in)      :: this
        type(c_ptr), intent(inout)           :: df
        logical, optional, intent(in)        :: rescale_forward
        logical, optional, intent(in)        :: synchronize 

        logical   :: rescale_forward_local, synchronize_local
        integer   :: df_size 
        real(r64) :: norm_cnt
        complex(c_double_complex), pointer :: fortran_df_double(:)
        complex(c_float_complex),  pointer :: fortran_df_single(:)

        rescale_forward_local = .true.
        synchronize_local     = .true.

        if (present(rescale_forward)) rescale_forward_local = rescale_forward
        if (present(synchronize)) synchronize_local = synchronize

        df_size  = product(this%dims)
        norm_cnt = 1.0_r64 / df_size

        if (this%is_double) then
            call c_f_pointer(df, fortran_df_double, [df_size])
            call fftw_execute_dft(this%plan, fortran_df_double, fortran_df_double)

            if (this%fft_sign == FFTW_FORWARD .and. rescale_forward_local) then
                fortran_df_double(:) = norm_cnt * fortran_df_double(:)
            end if

            nullify(fortran_df_double)
        else
            call c_f_pointer(df, fortran_df_single, [df_size])
            call fftwf_execute_dft(this%plan, fortran_df_single, fortran_df_single)

            if (this%fft_sign == FFTW_FORWARD .and. rescale_forward_local) then
                fortran_df_single(:) = norm_cnt * fortran_df_single(:)
            end if

            nullify(fortran_df_single)
        end if

    end subroutine execute

    !> Cleans up the fft_device_t
    !> @param[in,out] this - the fft_device_t to clean up
    subroutine delete(this)
        class(fft_device_t), intent(inout) :: this

        if (this%is_double) then
            call fftw_destroy_plan(this%plan)
        else
            call fftwf_destroy_plan(this%plan)
        end if

        ! Clean the pointer to device world handler
        nullify(this%world)

    end subroutine delete

end module m_fft_device
