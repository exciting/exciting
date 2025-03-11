! Copyright (C) 2024 exciting developers
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
!
!> This module contains the global object that
!> controls communication between device-target
module mod_device_offload
    use iso_c_binding,    only: c_int
    use iso_fortran_env,  only: i32 => int32
    use m_device_world_t, only: device_world_t

    implicit none


    private
    public :: init_device_world, &
              finish_device_world, &
              device_world

    !> This object controls the communication between
    !> the host and the device
    type(device_world_t) :: device_world

contains

    !> This subroutine inits the device world and other related quantities
    !> It associates all processes to a device
    subroutine init_device_world(mpiworld)
        !> The MPI communicator
        integer(i32), intent(in) :: mpiworld
        call device_world%init(int(mpiworld, kind=c_int))
    end subroutine init_device_world

    !> This subroutine finishes the device world
    subroutine finish_device_world()
        call device_world%synchronize()
        call device_world%finish()
    end subroutine finish_device_world

end module mod_device_offload
