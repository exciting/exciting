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
    contains
        procedure, public :: init, finish, is_queue_set, get_queue, syncronize, get_device, get_num_teams
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

        integer :: num_teams
        integer :: nprocs, myrank, ierr
        integer :: device
        integer, allocatable :: host_ids(:), device_ids(:)
        logical, allocatable :: local(:)
        integer :: i
        integer(c_int) :: cerror

        ! Get MPI information
        call mpi_comm_size(world, nprocs, ierr)
        call mpi_comm_rank(world, myrank, ierr)
        ! Allocate arrays holding rank dependent information
        allocate(local(nprocs), host_ids(nprocs), device_ids(nprocs))

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Init MAGMA
        call magma_init()
#endif  
        ! The omp_get_initial_device does not provide unique tag for 
        ! hosts in the MPI framework, i.e. two ranks in different
        ! nodes can have the same id. Nevertheless, that is internally
        ! used by the register as it is the one required by OpenMP routines
        this%host = get_host_id()
        ! Get a list of all procs
        call mpi_allgather(this%host, 1, MPI_INTEGER, host_ids, 1, MPI_INTEGER, world, ierr)

        ! Determine which processors are on this node so can control one of its GPUs
        local(1:nprocs) = host_ids(:) == this%host
        ! This construction ensures that for systems with embeded cards
        ! those are discarded
        device = omp_get_default_device()
        device_ids(1:nprocs) = -42
        do i = 1, nprocs
            if (local(i)) then
                device_ids(i) = device
                device = device + 1
            end if
        end do

        ! Having multiple processes share devices is not recommended.
        ! Therefore, we fail if such is the case
        ! This construction ensures that for systems with embeded cards
        ! those are discarded
        this%ndevices = omp_get_num_devices() - omp_get_default_device()
        if (any(device_ids >= this%ndevices)) then
            error stop "Error(device_world_t%init): Having multiple processes share devices is not recommended."
        end if
        
        ! Set the device associated to the process
        ! This propagates to all GPU call from this process
        ! Except if they specifically change the GPU id
        this%device = device_ids(myrank+1)
        call omp_set_default_device(this%device)

#if defined(NVIDIAGPU) || defined(AMDGPU)
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
        !$omp target teams map(from: num_teams)
        num_teams   = omp_get_num_teams()
        !$omp end target teams
        this%num_teams   = num_teams
        this%num_threads = omp_get_max_threads()
    
        ! Init register
        call this%register%init()

        do i = 1, nprocs
            if (myrank == i-1) then
                write(*,*) 'GPU world information: rank (', myrank ,')'
                write(*,*) 'Host ', this%host
                write(*,*) 'Device id', this%device
                write(*,*) 'Teams', this%num_teams
                write(*,*) 'Threads', this%num_threads
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

end module m_device_world_t
