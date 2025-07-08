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
        !> Device queue for linear algebra
        type(c_ptr), allocatable, private  :: queue(:)
        !> The number of teams in the device
        integer, private :: num_teams
        !> The maximum number of threads per team
        integer, private :: num_threads
        !> Flag to indicate if CPU-only backend is used
        logical, private :: cpu_backend = .false.
    contains
        procedure, public :: init, finish, is_queue_set, get_queue, synchronize, get_device, get_num_teams, using_cpu_backend, &
                             get_num_threads, simd_size, get_stream
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
        integer :: fdevreport

        integer :: i

        character(len=*), parameter :: logname = "GPU.LOG"

        ! Get MPI information
        call mpi_comm_rank(world, myrank, ierr)
        call mpi_comm_size(world, nprocs, ierr)

        ! Get unique ID for the host
        this%host = get_host_id()

        if (myrank == 0) call execute_command_line("echo '# Device-Host report' > " // logname)

        ! Check if the GPU-CPU association is done through environment variables
        do i = 1, nprocs
            if (myrank == i-1) then
                open(newunit=fdevreport, file=logname, action='write', position='append')
#if defined(NVIDIAGPU)
                call get_environment_variable("CUDA_VISIBLE_DEVICES", value=visible_devices, status=set_in_environment)
                write(fdevreport,*) "GPU Init : Capturing CUDA_VISIBLE_DEVICES (rank =", myrank , ") :", trim(visible_devices) 
#endif
#if defined(AMDGPU)
                call get_environment_variable("ROCR_VISIBLE_DEVICES", value=visible_devices, status=set_in_environment)
                write(fdevreport,*) "GPU init : Capturing ROCR_VISIBLE_DEVICES (rank =", myrank , ") :", trim(visible_devices)
#endif
#if defined(INTELGPU)
                call get_environment_variable("OMP_DEFAULT_DEVICE", value=visible_devices, status=set_in_environment)
                write(fdevreport,*) "GPU init : Capturing OMP_DEFAULT_DEVICE (rank =", myrank , ") :", trim(visible_devices)
#endif
                close(fdevreport)
            end if
            call mpi_barrier(world, ierr)
        end do

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

            this%device    = 0 
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

            this%device    = 0 
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

        ! Allocating the queue
        allocate(this%queue(omp_get_max_threads()))

        ! For some reason out of my grasp, GNU fails for allocate(this%queue(omp_get_max_threads(), source=c_null_ptr))
        do i = 1, omp_get_max_threads()
            this%queue(i) = c_null_ptr
        end do

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Set default device in MAGMA
        call magma_set_device(this%device)
        ! Init MAGMA
        call magma_init()
        ! Init the MAGMA queue, notice each OpenMP thread will have a queue, thus allowing for concurrency 
        ! of small linear algebra kernels
        !$omp parallel shared(this)
        call magma_queue_create(this%device, this%queue(omp_get_thread_num()+1))
        !$omp end parallel
#endif

        ! For AMD cards we need to make a global init
#if defined(AMDGPU)
        cerror = rocfft_setup()
        if (cerror /= rocfft_status_success) then
            error stop "Error(device_world_t%init): AMD FFT library (rocFFT) initialization failed."
        end if
#endif
        ! Get the number of teams and threads in the device regions
        ! The 50000 is only to put some virtual pressure on the GPU. It should be 
        ! always larger than the "treads" of a GPU "execution unit" (i.e. 
        ! the bunch of threads that share fast-memory).
        !$omp target map(from: num_teams, num_threads) 
        !$omp teams distribute parallel do simd lastprivate(num_teams, num_threads)
        do i = 1, 50000
            num_teams   = omp_get_num_teams()
            num_threads = omp_get_num_threads()
        end do
        !$omp end teams distribute parallel do simd
        !$omp end target

        this%num_teams   = num_teams
        this%num_threads = num_threads

        ! Print info
        call mpi_comm_size(world, nprocs, ierr)

        do i = 1, nprocs
            if (myrank == i-1) then
                open(newunit=fdevreport, file=logname, action='write', position='append')
                write(fdevreport,*) 'GPU world information: rank (', myrank ,')'
                write(fdevreport,*) 'Host ', this%host
                write(fdevreport,*) 'Device id ', this%device, 'of', this%ndevices
                write(fdevreport,*) 'Teams ', this%num_teams
                write(fdevreport,*) 'Threads ', this%num_threads
                write(fdevreport,*) 'Automatic host-device association ', set_in_environment == 1
                close(fdevreport)
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
        call this%synchronize()

        ! Final call for FFT libraries
#if defined(AMDGPU)
        cerror = rocfft_cleanup()
        if (cerror /= rocfft_status_success) then
            error stop "Error(device_world_t%finish): AMD FFT library (rocFFT) cleanup failed."
        end if
#endif

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Destroy queue
        !$omp parallel shared(this)
        call magma_queue_destroy(this%queue(omp_get_thread_num()+1))
        !$omp end parallel

        ! Destroy world
        call magma_finalize()
#endif
        
        ! Deallocate queue array
        deallocate(this%queue)

    end subroutine finish

    !> This function returns true if the queue is inited
    !> @param[in] this - the world to check if has inited queue
    function is_queue_set(this) result(answer)

        class(device_world_t), intent(in) :: this
        logical :: answer

        answer = C_associated(this%queue(omp_get_thread_num()+1))

    end function is_queue_set

    !> This is a getter for the queue
    !> @param[in] this - the world object from which the queue is retrieved
    function get_queue(this) result(queue)

        class(device_world_t), target, intent(in) :: this
        type(c_ptr) :: queue
        
        queue = this%queue(omp_get_thread_num()+1)

    end function get_queue

    !> This synchronizes the world 
    !> @param[in] this - the device which we want to sync with
    subroutine synchronize(this)
        class(device_world_t), target, intent(in) :: this

        integer(c_int) :: cerror

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! This synchronizes the MAGMA stream, which is also used for FFTs
        call magma_queue_sync(this%queue(omp_get_thread_num()+1))
#endif
        ! We also wait for all nonwait OpenMP pure kernels
        !$omp taskwait
    end subroutine synchronize

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

    !> Returns the underlying stream for control
    !> for Intel returns the pointer of the queue, which can be used for
    !> depend constructs
    type(c_ptr) function get_stream(this)
        class(device_world_t), intent(in) :: this
#if defined(NVIDIAGPU) 
        get_stream = magma_queue_get_cuda_stream(this%queue(omp_get_thread_num()+1))
#elif defined(AMDGPU)
        get_stream = magma_queue_get_hip_stream(this%queue(omp_get_thread_num()+1))
#else
        get_stream = this%queue(omp_get_thread_num()+1)
#endif
    end function get_stream

end module m_device_world_t
