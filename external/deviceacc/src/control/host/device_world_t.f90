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
!> This file contains things to handle host-device communication. CPU backend.

!> Module containing things to handle host-device communication
module m_device_world_t
    
    use iso_c_binding

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
        integer, private :: num_teams = 1
        !> The maximum number of threads per team
        integer, private :: num_threads = 1 
        !> Flag to indicate if CPU-only backend is used
        logical, private :: cpu_backend = .true.
    contains
        procedure, public :: init, finish, is_queue_set, get_queue, synchronize, get_device, get_num_teams, using_cpu_backend, &
                             get_num_threads, simd_size, get_stream
    end type device_world_t

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
    end subroutine init

    !> This subroutine finishes the device handler
    !> @param[in] this - the device handler to finish
    subroutine finish(this)
        class(device_world_t), intent(inout) :: this
    end subroutine finish

    !> This function returns true if the queue is inited
    !> @param[in] this - the world to check if has inited queue
    function is_queue_set(this) result(answer)
        class(device_world_t), intent(in) :: this
        logical :: answer
        answer = .false.
    end function is_queue_set

    !> This is a getter for the queue
    !> @param[in] this - the world object from which the queue is retrieved
    function get_queue(this) result(queue)
        class(device_world_t), target, intent(in) :: this
        type(c_ptr) :: queue
        queue = c_null_ptr
    end function get_queue

    !> This synchronizes the world (in a very agresive way)
    !> @param[in] this - the device which we want to sync with
    subroutine synchronize(this)
        class(device_world_t), target, intent(in) :: this
        !$omp barrier
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
        simd_size = 1
    end function simd_size

    !> Returns the underlying stream, or object for depend clauses
    type(c_ptr) function get_stream(this)
        class(device_world_t), intent(in) :: this
        get_stream = c_null_ptr
    end function get_stream

end module m_device_world_t
