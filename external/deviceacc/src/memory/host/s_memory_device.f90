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
!> Contains the CPU backend to allocations
submodule(m_memory_device) sm_memory_device

    implicit none

contains

    module subroutine allocate_device_memory(memory, memsize, device_id)
        type(c_ptr), intent(out) :: memory 
        integer(c_size_t), intent(in) :: memsize
        integer, intent(in) :: device_id

        memory = c_malloc(memsize)

    end subroutine allocate_device_memory

    module subroutine deallocate_device_memory(memory, device_id)
        type(c_ptr), intent(inout) :: memory 
        integer, intent(in) :: device_id

        call c_free(memory)

    end subroutine deallocate_device_memory

    module type(c_ptr) function get_device_pointer(host_data, device_id)
        type(*), target, intent(in) :: host_data(..)
        integer, intent(in) :: device_id
        get_device_pointer = c_loc(host_data)
    end function get_device_pointer

end submodule sm_memory_device
