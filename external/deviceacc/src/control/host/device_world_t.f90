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
    use omp_lib
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
        ! Init register
        call this%register%init()
    end subroutine init

    !> This subroutine finishes the device handler
    !> @param[in] this - the device handler to finish
    subroutine finish(this)
        class(device_world_t), intent(inout) :: this
        ! Clean register
        call this%register%finish()
    end subroutine finish

    !> This function returns true if the queue is inited
    !> @param[in] this - the world to check if has inited queue
    pure function is_queue_set(this) result(answer)
        class(device_world_t), intent(in) :: this
        logical :: answer
        answer = .false.
    end function is_queue_set

    !> This is a getter for the queue
    !> @param[in] this - the world object from which the queue is retrieved
    function get_queue(this) result(queue)
        class(device_world_t), target, intent(in) :: this
        type(C_ptr), pointer :: queue
        queue = c_null_ptr
    end function get_queue

    !> This syncronizes the world (in a very agresive way)
    !> @param[in] this - the device which we want to sync with
    subroutine syncronize(this)
        class(device_world_t), target, intent(in) :: this
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
