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
!> Supported vendors include NVIDIA, Intel, and AMD. CPU-only backend.
module m_fft_device
    
    use iso_c_binding,    only: c_f_pointer, c_int, c_ptr, c_float_complex, c_loc, &
                                c_double_complex, c_size_t, c_null_ptr 
    use iso_fortran_env,  only: i32=>int32, r32=>real32, r64=>real64
    use m_device_world_t, only: device_world_t

#if defined(AMDGPU)
    use hipfort_rocfft, only: rocfft_precision_double, rocfft_precision_single, &
                              rocfft_transform_type_complex_forward, rocfft_transform_type_complex_inverse, &
                              rocfft_placement_inplace, rocfft_plan_create, rocfft_status_success, &
                              rocfft_execution_info_create, rocfft_execution_info_set_stream, &
                              rocfft_execute, rocfft_execution_info_destroy, rocfft_plan_destroy
#endif 
#if defined(INTELGPU)
    use mkl_dft_type, only: dfti_descriptor, DFTI_DOUBLE, DFTI_SINGLE, DFTI_FORWARD_SCALE, DFTI_COMPLEX 
    use mkl_dfti_omp_offload, only: dfti_create_descriptor_highd, dfti_compute_forward_z_cpu, &
                                    dfti_compute_backward_z_cpu, dfti_compute_forward_c_cpu, &
                                    dfti_compute_backward_c_cpu, DftiFreeDescriptor, DftiSetValue
#endif

    implicit none

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
#if defined(INTELGPU)
        !> Handler for FFT
        type(dfti_descriptor), private, pointer :: descriptor
#endif 
        !> Actual FFT plan.
#if defined(NVIDIAGPU)
        integer(c_int), private :: plan 
#else 
        type(c_ptr), private :: plan = c_null_ptr
#endif
        !> Double precision
        logical, private :: is_double
        !> Device world
        type(device_world_t), pointer :: world
    contains
        procedure :: initialize
        procedure :: delete
        procedure :: execute
    end type fft_device_t

    
#if defined(NVIDIAGPU)
    ! NVIDIA cuFFT interface (this removes the need for CUDA Fortran but uses C CUDA api
    ! so no nvfortran is needed. Be aware that in cuFFT only up to 3D, trying higher dimensions will
    ! fail. 
interface

        !> Interface to CUDA FFT plan destructor C function
        !> @parameter[in,out] plan - the fft plan
        function cufftDestroy(plan) result(info) bind (c, name="cufftDestroy")
            import
            integer(c_int), value :: plan
            integer(c_int)        :: info
        end function cufftDestroy

        !> Interface to 1D CUDA FFT plan creator C function
        !> @parameter[out] plan - the fft plan
        !> @parameter[in]  nx   - first dimension size (C-order) 
        !> @parameter[in]  type - flag to determine what to do
        !> @parameter[in]  batch - the number of FFT to do
        !> @result         info - the info about the execution
        function cufftPlan1d(plan, nx, type, batch) result(info) bind (c, name="cufftPlan1d")
            import
            integer(c_int), value :: nx
            integer(c_int), value :: type
            integer(c_int), value :: batch
            integer(c_int), intent(out) :: plan
            integer(c_int) :: info
        end function cufftPlan1d

        
        !> Interface to 2D CUDA FFT plan creator C function
        !> @parameter[out] plan - the fft plan
        !> @parameter[in]  nx   - first dimension size (C-order) 
        !> @parameter[in]  ny   - 2nd dimension size (C-order)
        !> @parameter[in]  type - flag to determine what to do
        !> @result         info - the info about the execution
        function cufftPlan2d(plan, nx, ny, type) result(info) bind (c, name="cufftPlan2d")
            import
            integer(c_int), value :: nx
            integer(c_int), value :: ny
            integer(c_int), value :: type
            integer(c_int), intent(out) :: plan
            integer(c_int) :: info
        end function cufftPlan2d

        !> Interface to 3D CUDA FFT plan creator C function
        !> @parameter[out] plan - the fft plan
        !> @parameter[in]  nx   - first dimension size (C-order) 
        !> @parameter[in]  ny   - 2nd dimension size (C-order)
        !> @parameter[in]  nz   - 3rd dimension size (C-order)
        !> @parameter[in]  type - flag to determine what to do
        !> @result         info - the info about the execution 
        function cufftPlan3d(plan, nx, ny, nz, type) result(info) bind (c, name="cufftPlan3d")
            import
            integer(c_int), value :: nx
            integer(c_int), value :: ny
            integer(c_int), value :: nz
            integer(c_int), value :: type
            integer(c_int), intent(out) :: plan
            integer(c_int) :: info
        end function cufftPlan3d

        !> Executes a CUDA FFT plan
        !> @parameter[in] plan - the fft plan
        !> @parameter[in] in - pointer to the input data
        !> @parameter[out] out - pointer to where the results are to be stored
        !> @parameter[in] direction - the direction of the transform
        function cufftXtExec(plan, in, out, direction) result(info) bind(c, name="cufftXtExec")
            import
            integer(c_int), value :: plan
            type(c_ptr), value    :: in 
            type(c_ptr), value    :: out 
            integer(c_int), value :: direction
            integer(c_int) :: info
        end function  cufftXtExec

        !> Associates a plan to a stream
        function cufftSetStream(plan, stream) result(info) bind(c, name="cufftSetStream")
            import
            integer(c_int), value :: plan
            type(c_ptr), value    :: stream
            integer(c_int) :: info
        end function cufftSetStream

end interface
! END if defined(NVIDIAGPU)
#endif

contains 

    !> Initializes the FFT type using device accelerated procedures
    !> @param[out] this    - the fft_device_t to initialize
    !> @param[in]  dims    - the shape of the problem to solve
    !> @param[in]  fft_sign    - the direction of the Fourier transform (forward: -1 // backward: 1) 
    !> @param[in,out]  df      - the data or that of similar shape over which the FFT is computed (device pointer)
    !> @param[in]  is_double  - are we dealing with double precision, be aware as we use c_ptr the type information is lost
    !> @param[in,out] world - device handler 
    subroutine initialize(this, dims, fft_sign, df, is_double, world)
        
        class(fft_device_t), intent(out)            :: this
        integer, intent(in)                         :: dims(:)
        integer, intent(in)                         :: fft_sign
        type(c_ptr), intent(inout)                  :: df
        logical, intent(in)                         :: is_double 
        type(device_world_t), target, intent(inout) :: world 

        integer(c_int) :: error
        integer(c_int) :: fft_precision, fft_task
        integer(c_size_t), allocatable, target :: amd_dims(:)

        this%rank   = int(size(dims), kind = c_int)
#if defined(NVIDIA) || defined(AMDGPU)
        this%dims   = int(dims(this%rank:1:-1), kind = c_int) ! Used libraries expect C-order
#endif
#if defined(INTELGPU)
        this%dims = int(dims, kind = c_int) 
#endif
        this%fft_sign   = fft_sign
        this%is_double = is_double
        this%world  => world

        ! Some sanity check
        if (this%fft_sign /= -1 .and. this%fft_sign /= 1) then
            error stop "fft_device_t%initialize: The only accepted fft_signs are -1 (forward) and 1 (backward)."
        end if

#if defined(NVIDIAGPU)
        ! Get precision, the flag is an integer with values written as hexadecimal numbers
        if (this%is_double) then
            fft_precision = int(z'69')
        else
            fft_precision = int(z'29')
        end if

        select case(this%rank)
        case(1)
            error = cufftPlan1d(this%plan, this%dims(1), fft_precision, 1_c_int)
        case(2)
            error = cufftPlan2d(this%plan, this%dims(1), this%dims(2), fft_precision)
        case(3)
            error = cufftPlan3d(this%plan, this%dims(1), this%dims(2), this%dims(3), fft_precision)
        case default
            error stop "fft_device_t%initialize: cuFFT does only support 1, 2 and 3D FFT."
        end select

        if (error /= 0) then
            error stop "Error(fft_device_t%initialize): CUDA FFT failed to create a plan"
        end if

        ! Now associate the plan with the stream
        error = cufftSetStream(this%plan, this%world%get_stream())
        if (error /= 0) then
            error stop "Error(fft_device_t%initialize): CUDA FFT failed to associate stream to the plan"
        end if
! END if defined(NVIDIAGPU)
#endif
#if defined(AMDGPU)
        ! AMD needs to know the precision in each plan
        if (this%is_double) then
            fft_precision = rocfft_precision_double
        else
            fft_precision = rocfft_precision_single
        end if
        ! AMD needs to know what to do
        if (this%fft_sign == -1) then
            fft_task = rocfft_transform_type_complex_forward
        else
            fft_task = rocfft_transform_type_complex_inverse
        end if
        
        amd_dims = int(this%dims, c_size_t)

        ! Call to rocFFT function to create the plan
        error = rocfft_plan_create(this%plan, &
                                   rocfft_placement_inplace, &
                                   fft_task, &
                                   fft_precision, &
                                   int(this%rank, c_size_t), &
                                   c_loc(amd_dims), &
                                   1_c_size_t, &
                                   c_null_ptr)

        if (error /= rocfft_status_success) then
            error stop "Error(fft_device_t%initialize): rocfft_plan_create failed"
        end if
! END if defined(AMDGPU)
#endif

#if defined(INTELGPU)
        ! Using the native MKL DFTI
        if (this%is_double) then
            fft_precision = DFTI_DOUBLE
        else
            fft_precision = DFTI_SINGLE
        end if
        ! Creating a handler
        ! Notice that the handler does not need the array in this case
        error = dfti_create_descriptor_highd(this%descriptor, fft_precision, DFTI_COMPLEX, this%rank, this%dims)
        ! Checking 
        if (error /= 0) then
            error stop "Error(fft_device_t%initialize): dfti_create_descriptor_highd failed"
        end if
#endif

    end subroutine initialize

    !> Executes the plan
    !> @param[in,out]      this - the fft_device_t for which to execute the plan
    !> @param[in,out]  df   - the data over which to perform the plan (device ptr). On exit contains the result.
    !> @param[in]      rescale_forward - rescale the FFT in case of forward FFT
    !> @param[in]      synchronize - force syncronization after the execution of the plan
    subroutine execute(this, df, rescale_forward, synchronize)

        class(fft_device_t), intent(in)      :: this
        type(c_ptr), intent(inout)           :: df
        logical, optional, intent(in)        :: rescale_forward
        logical, optional, intent(in)        :: synchronize 

        logical   :: rescale_forward_local, synchronize_local
        integer   :: i, df_size 
        real(r64) :: norm_cnt
        complex(c_double_complex), pointer :: fortran_df_double(:)
        complex(c_float_complex),  pointer :: fortran_df_single(:)
        integer(c_int) :: error
        complex(r64), target, allocatable :: tin(:), tout(:)
        type(c_ptr) :: execution_info
        
        rescale_forward_local = .true.
        synchronize_local      = .true.

        if (present(rescale_forward)) rescale_forward_local = rescale_forward
        if (present(synchronize)) synchronize_local = synchronize

        ! Execute the plan
#if defined(AMDGPU)
        ! We will associate the execution to the current stream
        error = rocfft_execution_info_create(execution_info)
        if (error /= rocfft_status_success) then
            error stop "Error(fft_device_t%execute): rocfft_execution_info_create failed"
        end if
 
        error = rocfft_execution_info_set_stream(execution_info, this%world%get_stream())
        if (error /= rocfft_status_success) then
            error stop "Error(fft_device_t%execute): rocfft_execution_info_set_stream failed"
        end if

        error = rocfft_execute(this%plan, df, c_null_ptr, execution_info)
        if (error /= rocfft_status_success) then
            error stop "Error(fft_device_t%execute): rocfft_execute failed"
        end if

        error = rocfft_execution_info_destroy(execution_info)
        if (error /= rocfft_status_success) then
            error stop "Error(fft_device_t%execute): rocfft_execution_info_destroy failed"
        end if

#endif 

#if defined(NVIDIAGPU)
        ! Notice that we do it inplace
        error = cufftXtExec(this%plan, df, df, int(this%fft_sign, kind = c_int))
        if (error /= 0) then
            error stop "Error(fft_device_t%execute): cufftXtExec"
        end if
#endif
        ! In the case of NVIDIA and AMD we use an offloaded 
        ! do loop to rescale within the device the quantities if asked 
#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Syncronize the host and the device 
        if (synchronize_local) then 
            call this%world%synchronize()
        end if

        ! Rescale, be aware that the rescale is done in the device data 
        if (this%fft_sign == -1 .and. rescale_forward_local) then
            
            ! We need an explicit check regarding if FFT is done 
            if (.not. synchronize_local) then 
                call this%world%synchronize()
            end if

            df_size  = product(this%dims)
            ! Either for device or host multiplication is faster than division
            norm_cnt = 1.0_r64 / df_size
            
            ! The pointers point to device memory
            ! so the data will remain on the device during execution of the device-offloaded region
            ! That means no transfer happens
            if (this%is_double) then
                call c_f_pointer(df, fortran_df_double, [df_size])
                !$omp target teams distribute parallel do simd private(i) default(shared)
                do i = 1, df_size
                    fortran_df_double(i) = fortran_df_double(i) * norm_cnt
                end do
                !$omp end target teams distribute parallel do simd
                nullify(fortran_df_double) 
            else 
                call c_f_pointer(df, fortran_df_single, [df_size])
                !$omp target teams distribute parallel do simd private(i) default(shared)
                do i = 1, df_size
                    fortran_df_single(i) = fortran_df_single(i) * norm_cnt
                end do
                !$omp end target teams distribute parallel do simd
                nullify(fortran_df_single) 
            end if
        end if
! END if defined(NVIDIAGPU) || defined(AMDGPU)
#endif
        
#if defined(INTELGPU)

        ! Note that current dispatch in interfaces blocks cause in ifx v 2024.0.2
        ! the compilation to break. Therefore we do not use the interface but the specific calls
        ! TODO: revise in future versions, so we can remove the specific calls and use the interface

        df_size  = product(this%dims)

        ! Select the acording function depending on the direction of the FFT
        if (this%fft_sign == -1) then
            
            ! Either for device or host multiplication is faster than division
            norm_cnt = 1.0_r64 / df_size

            ! MKL DFTI rescaling during the execution
            if (this%fft_sign == -1 .and. rescale_forward_local) then
                error = DftiSetValue(this%descriptor, DFTI_FORWARD_SCALE, norm_cnt)
                error stop "Error(fft_device_t%execute): DftiSetValue failed"
            end if

            ! Depending on the kind of data fill one or the other pointer
            if (this%is_double) then
                call c_f_pointer(df, fortran_df_double, [df_size])
                !$omp dispatch is_device_ptr(fortran_df_double)
                error = dfti_compute_forward_z_cpu(this%descriptor, fortran_df_double)
                nullify(fortran_df_double)
            else 
                call c_f_pointer(df, fortran_df_single, [df_size])
                !$omp dispatch is_device_ptr(fortran_df_single)
                error = dfti_compute_forward_c_cpu(this%descriptor, fortran_df_single)
                nullify(fortran_df_single)
            end if

            if (error /= 0) then
                error stop "Error(fft_device_t%execute): DftiComputeForward failed"
            end if

        else            
            ! Depending on the kind of data fill one or the other pointer 
            if (this%is_double) then
                call c_f_pointer(df, fortran_df_double, [df_size])
                !$omp dispatch is_device_ptr(fortran_df_double)
                error = dfti_compute_backward_z_cpu(this%descriptor, fortran_df_double)
                nullify(fortran_df_double)
            else 
                call c_f_pointer(df, fortran_df_single, [df_size])
                !$omp dispatch is_device_ptr(fortran_df_single)
                error = dfti_compute_backward_c_cpu(this%descriptor, fortran_df_single)
                nullify(fortran_df_single)
            end if
            
            if (error /= 0) then
                error stop "Error(fft_device_t%execute): DftiComputeBackward failed"
            end if

        end if

        ! Syncronize the host and the device 
        if (synchronize_local) call this%world%synchronize()
! END if defined(INTELGPU)
#endif

    end subroutine execute

    !> Cleans up the fft_device_t
    !> @param[in,out] this - the fft_device_t to clean up
    subroutine delete(this)
        class(fft_device_t), intent(inout) :: this
        
        integer(c_int) :: error 

#if defined(AMDGPU)
        error = rocfft_plan_destroy(this%plan)
        if (error /= rocfft_status_success) then
            error stop "Error(fft_device_t%delete): rocfft_plan_destroy failed"
        end if
#endif
#if defined(NVIDIAGPU)
        error = cufftDestroy(this%plan)
        if (error /= 0) then
            error stop "Error(fft_device_t%delete): cufftDestroy failed"
        end if
#endif
#if defined(INTELGPU)
        ! Clean the descriptor
        error = DftiFreeDescriptor(this%descriptor)
#endif
        ! Clean the pointer to device world handler
        nullify(this%world)

    end subroutine delete

end module m_fft_device
