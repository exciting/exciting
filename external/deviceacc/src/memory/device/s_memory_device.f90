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
!> Contains the device backend to allocations
submodule(m_memory_device) sm_memory_device

    use omp_lib, only: omp_target_alloc, omp_target_free, omp_get_mapped_ptr

    implicit none

contains

    module subroutine allocate_device_memory(memory, memsize, device_id)
        type(c_ptr), intent(out) :: memory 
        integer(c_size_t), intent(in) :: memsize
        integer, intent(in) :: device_id

#if !defined(_USM_)
        memory = omp_target_alloc(memsize, device_id)
#else
        memory = c_malloc(memsize)
#endif

        if (.not. c_associated(memory)) error stop "Error(allocate_device_memory) : returned null pointer"

    end subroutine allocate_device_memory

    module subroutine deallocate_device_memory(memory, device_id)
        type(c_ptr), intent(inout) :: memory 
        integer, intent(in) :: device_id

        if (.not. c_associated(memory)) error stop "Error(deallocate_device_memory) : cannot free a null pointer"

#if !defined(_USM_)
        call omp_target_free(memory, device_id)
#else
        call c_free(memory)
#endif

    end subroutine deallocate_device_memory

    module type(c_ptr) function get_device_pointer(host_data, device_id)
        type(*), target, intent(in) :: host_data(..)
        integer, intent(in) :: device_id

#if !defined(_USM_)
        get_device_pointer = omp_get_mapped_ptr(c_loc(host_data), device_id)
#else
        get_device_pointer = c_loc(host_data)
#endif
        
        if (.not. c_associated(get_device_pointer)) error stop "Error(get_device_pointer) : returned null pointer"
    
    end function get_device_pointer

end submodule sm_memory_device
