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
!> This file contains things to handle host-device communication. GPU backend.

!> Module containing things to handle host-device communication
#if defined(INTELGPU)
include "mkl_omp_offload.f90"
#endif
module m_device_world_t
    
    use iso_c_binding
    use omp_lib
#if defined(NVIDIAGPU)
    use magma2
#endif
#if defined(AMDGPU)
    use magma2
    use hipfort_rocfft
#endif
    use m_device_host_register_fortran, only: device_host_register
    use mpi

    implicit none

    private
    public  device_world_t
    
    !> Type to handle the devices
    type device_world_t
        !> Device id taking care of control
        integer, private :: device   = -1
        !> Host id 
        integer(c_int), private :: host = -1
        !> Number of devices
        integer, private :: ndevices = -1
        !> Device queue for MAGMA
        type(c_ptr), private  :: queue = c_null_ptr
        !> The number of teams in the device
        integer, private :: num_teams
        !> The maximum number of threads per team
        integer, private :: num_threads
        !> Device host register
        type(device_host_register) :: register
        !> Flag to indicate if CPU-only backend is used
        logical, private :: cpu_backend = .false.
    contains
        procedure, public :: init, finish, is_queue_set, get_queue, syncronize, get_device, get_num_teams, using_cpu_backend, &
                             get_num_threads, simd_size, get_linalg_stream
    end type device_world_t

interface
    !> Gets the unique identifier of the current host
    !> @result get_host_id - the host unique id
    function get_host_id() bind(C, name="gethostid")
      import 
      integer (c_int) :: get_host_id
    end function get_host_id
end interface

#if defined(AMDGPU)
interface
    function hipDeviceSynchronize() bind(c, name="hipDeviceSynchronize")
        use iso_c_binding
        integer(c_int) :: hipDeviceSynchronize
    end function hipDeviceSynchronize

#if defined(AMD_CAN_SET_VALID_DEVICES)
    function hipSetValidDevices(device_arr, len) bind(c, name="hipSetValidDevices")
        use iso_c_binding
        integer(c_int) :: hipSetValidDevices  ! Return type
        integer(c_int), intent(in) :: device_arr(len)  ! Array of device IDs
        integer(c_int), value :: len  ! Length of the array
    end function hipSetValidDevices
#endif

end interface
#endif

#if defined(NVIDIAGPU) 
interface 
    function cudaSetValidDevices(device_arr, len) bind(c, name="cudaSetValidDevices")
        use iso_c_binding
        integer(c_int) :: cudaSetValidDevices  ! Return type
        integer(c_int), intent(in) :: device_arr(len)  ! Array of device IDs
        integer(c_int), value :: len  ! Length of the array
    end function cudaSetValidDevices
end interface 
#endif

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                   DEVICE WORLD                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> This subroutine inits the device world and other related quantities
    !> @param[in] this    - the GPU magma world to initialize
    !> @param[in] world   - the MPI communicator
    subroutine init(this, world)
            
        class(device_world_t), target, intent(inout) :: this
        integer(C_int), intent(in) :: world

        integer :: num_teams, num_threads
        integer :: nprocs, myrank, ierr
        integer(c_int)         :: device_node_id(1)
        integer(c_int)         :: cerror
        integer(c_int)         :: host_world
        character(len=256)     :: visible_devices
        integer :: set_in_environment

        integer :: i

        ! Get MPI information
        call mpi_comm_rank(world, myrank, ierr)

        ! Get unique ID for the host
        this%host = get_host_id()

        ! Check if the GPU-CPU association is done through environment variables
#if defined(NVIDIAGPU)
        call get_environment_variable("CUDA_VISIBLE_DEVICES", value=visible_devices, status=set_in_environment)
        write(*,*) "GPU Init : Capturing CUDA_VISIBLE_DEVICES (rank =", myrank , ") :", trim(visible_devices) 
#endif
#if defined(AMDGPU)
        call get_environment_variable("ROCR_VISIBLE_DEVICES", value=visible_devices, status=set_in_environment)
        write(*,*) "GPU init : Capturing ROCR_VISIBLE_DEVICES (rank =", myrank , ") :", trim(visible_devices)
#endif
#if defined(INTELGPU)
        call get_environment_variable("OMP_DEFAULT_DEVICE", value=visible_devices, status=set_in_environment)
        write(*,*) "GPU init : Capturing OMP_DEFAULT_DEVICE (rank =", myrank , ") :", trim(visible_devices)
#endif

        ! If done by the user, the automatic association is ignored
        ! scan checks that visible_devices assingns only a host per device.
        ! TODO: In the future end.
        if (set_in_environment == 1 .or. scan(visible_devices, ",") /= 0) then
            ! Create a communicator between processes sharing a same physical node
            call mpi_comm_split_type(world, mpi_comm_type_shared, 0, mpi_info_null, host_world, ierr)
            ! Get the rank within the node
            call mpi_comm_rank(host_world, device_node_id(1), ierr)
            ! Remove embeded cards
            device_node_id(1) = device_node_id(1) + omp_get_default_device()

            ! Get number of processes within the node
            call mpi_comm_size(host_world, nprocs, ierr)

            ! Having multiple processes share devices is not recommended.
            ! Therefore, we fail in that case
            if ( nprocs > omp_get_num_devices() - omp_get_default_device() ) then
                error stop "Error(device_world_t%init): Having multiple processes share devices is not recommended."
            end if

#if defined(NVIDIAGPU)
            
            cerror = cudaSetValidDevices(device_node_id, 1_c_int)
            if (cerror /= 0) then
                error stop "Error(device_world_t%init): failed cudaSetValidDevices with "
            end if 

            if (omp_get_num_devices() /= 1) then
                error stop "Error(device_world_t%init): the number of devices has not properly limited"    
            end if

            this%device    = 0 !device_node_id(1)
            this%ndevices  = 1
#endif

#if defined(AMDGPU) && defined(AMD_CAN_SET_VALID_DEVICES)
            
            ! Note that this fails in LUMI currently but works in NVIDIA
            ! Ask AMD what would that mean (so if this would correspond to a runtime modification 
            ! of ROCR_VISIBLE_DEVICES)
            cerror = hipSetValidDevices(device_node_id, 1_c_int)
            if (cerror == 801) then
                error stop "Error(device_world_t%init): hipSetValidDevices unsupported."
            end if

            if (cerror /= 0) then
                error stop "Error(device_world_t%init): failed hipSetValidDevices"
            end if

            if (omp_get_num_devices() /= 1) then
                error stop "Error(device_world_t%init): the number of devices has not properly limited"
            end if

            this%device    = 0 !device_node_id(1)
            this%ndevices  = 1

#endif

#if ( defined(AMDGPU) && !defined(AMD_CAN_SET_VALID_DEVICES) ) || defined(INTELGPU)
            ! Restrict what the GPU world can see
            ! Set the device associated to the process
            this%device    = device_node_id(1)
            this%ndevices  = omp_get_num_devices()
#endif
        else
            this%device    = omp_get_default_device()
            this%ndevices  = omp_get_num_devices()
        end if

        ! Set OpenMP to use the selected device
        call omp_set_default_device(this%device)

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Set default device in MAGMA
        call magma_set_device(this%device)
        ! Init MAGMA
        call magma_init()
        ! Init the MAGMA queue
        call magma_queue_create(this%device, this%queue)
#endif

        ! For AMD cards we need to make a global init
#if defined(AMDGPU)
        cerror = rocfft_setup()
        if (cerror /= rocfft_status_success) then
            error stop "Error(device_world_t%init): AMD FFT library (rocFFT) initialization failed."
        end if
#endif
        ! Get the number of teams and threads in the device regions
#if defined(NVIDIAGPU) || defined(AMDGPU)  
        !$omp target map(from: num_teams, num_threads) 
        !$omp teams distribute parallel do
        do i = 1, 1
            num_teams   = omp_get_num_teams()
            num_threads = omp_get_num_threads()
        end do
        !$omp end teams distribute parallel do
        !$omp end target 
#endif
#if defined(INTELGPU)
        !$omp target map(from: num_teams, num_threads) 
        !$omp teams distribute parallel do
        do i = 1, 1
            num_teams   = omp_get_max_teams()
            num_threads = omp_get_max_threads()
        end do
        !$omp end teams distribute parallel do
        !$omp end target 
#endif
        this%num_teams   = num_teams
        this%num_threads = num_threads
        
        ! Init register
        call this%register%init()

        ! Print info
        call mpi_comm_size(world, nprocs, ierr)

        do i = 1, nprocs
            if (myrank == i-1) then
                write(*,*) 'GPU world information: rank (', myrank ,')'
                write(*,*) 'Host ', this%host
                write(*,*) 'Device id ', this%device, 'of', this%ndevices
                write(*,*) 'Teams ', this%num_teams
                write(*,*) 'Threads ', this%num_threads
                write(*,*) 'Automatic host-device association ', set_in_environment == 1
            end if
            call mpi_barrier(world, ierr)
        end do

    end subroutine init

    !> This subroutine finishes the device handler
    !> @param[in] this - the device handler to finish
    subroutine finish(this)
        
        class(device_world_t), intent(inout) :: this
        integer(c_int) :: cerror

        ! Make a final sync call only in case
        call this%syncronize()

        ! Final call for FFT libraries
#if defined(AMDGPU)
        cerror = rocfft_cleanup()
        if (cerror /= rocfft_status_success) then
            error stop "Error(device_world_t%finish): AMD FFT library (rocFFT) cleanup failed."
        end if
#endif

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Destroy queue
        call magma_queue_destroy(this%queue)
        
        ! Destroy world
        call magma_finalize()
#endif
        ! Clean register
        call this%register%finish()
    
    end subroutine finish

    !> This function returns true if the queue is inited
    !> @param[in] this - the world to check if has inited queue
    pure function is_queue_set(this) result(answer)

        class(device_world_t), intent(in) :: this
        logical :: answer

        answer = C_associated(this%queue)

    end function is_queue_set

    !> This is a getter for the queue
    !> @param[in] this - the world object from which the queue is retrieved
    function get_queue(this) result(queue)

        class(device_world_t), target, intent(in) :: this
        type(C_ptr), pointer :: queue
        
        queue => this%queue 

    end function get_queue

    !> This syncronizes the world (in a very agresive way)
    !> @param[in] this - the device which we want to sync with
    subroutine syncronize(this)
        class(device_world_t), target, intent(in) :: this

        integer(c_int) :: cerror

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_queue_sync(this%queue)
#endif
#if defined(AMDGPU)
        cerror = hipDeviceSynchronize()
        if (cerror /= 0_c_int) then
            error stop "Error(device_world_t%syncronize): error calling hipDeviceSynchronize"
        end if
#endif
        !$omp barrier
    end subroutine syncronize

    !> This provides device id
    !> @param[in] this - return the associated device id
    pure function get_device(this) result(device)
        class(device_world_t), intent(in) :: this
        integer :: device
        device = this%device
    end function get_device
    
    !> This provides the number of teams
    !> @param[in] this - return the number of teams of the device
    pure function get_num_teams(this) result(num_teams)
        class(device_world_t), intent(in) :: this
        integer :: num_teams
        num_teams = this%num_teams
    end function get_num_teams

    !> This provides the number of threads
    !> @param[in] this - return the number of teams of the device
    pure function get_num_threads(this) result(num_threads)
        class(device_world_t), intent(in) :: this
        integer :: num_threads
        num_threads = this%num_threads
    end function get_num_threads
    
    !> Returns .true. if using the CPU backend
    pure logical function using_cpu_backend(this)
        class(device_world_t), intent(in) :: this
        using_cpu_backend = this%cpu_backend
    end function using_cpu_backend

    !> Returns the size for SIMD in the device
    pure integer function simd_size(this)
        class(device_world_t), intent(in) :: this
#if defined(NVIDIAGPU)
        simd_size = 32
#endif
#if defined(AMDGPU)
        simd_size = 64
#endif
    end function simd_size

    !> Returns the underlying stream that handles linear algebra
    !> for Intel returns nothing
    type(c_ptr) function get_linalg_stream(this)
        class(device_world_t), intent(in) :: this
#if defined(NVIDIAGPU) 
        get_linalg_stream = magma_queue_get_cuda_stream(this%queue)
#elif defined(AMDGPU)
        get_linalg_stream = magma_queue_get_hip_stream(this%queue)
#else
        get_linalg_stream = c_null_ptr
#endif
    end function get_linalg_stream

end module m_device_world_t
