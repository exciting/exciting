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
!> This file contains wrappers to allocate/deallocate memory in the device using
!> c_ptr. These are the common elements, specific backends are compiled using submodules

! Classic ifort bug for versions 2021.4.x
#ifdef IFORT_2021_PURE_ASSUMED_RANK_BUG
#define IFORT_BUG_SAFE_PURE
#define IFORT_BUG_POLYCLASS_ASSUMED_RANK type(*)
#else
#define IFORT_BUG_SAFE_PURE pure
#define IFORT_BUG_POLYCLASS_ASSUMED_RANK class(*)
#endif

module m_memory_device

    use iso_c_binding, only: c_ptr, c_size_t, c_sizeof, c_float, & 
                             c_double, c_int, c_long, c_short, c_loc, &
                             c_null_ptr, c_associated, c_intptr_t

    implicit none

    private
    public :: allocate_device_memory, &
              deallocate_device_memory, &
              get_device_pointer, &
              bytes_size, &
              bytes_short, &
              bytes_int, &
              bytes_long_int, &
              bytes_single_real, &
              bytes_double_real, &
              bytes_single_complex, &
              bytes_double_complex, &
              get_pointer_alignment, &
              generate_batched_array

    !> The bit size of a bytes
    integer(c_size_t), parameter :: bit_size = storage_size('1', kind=c_size_t)
    !> The size in bytes of a short integer 
    integer(c_size_t), parameter :: bytes_short           = c_sizeof(1_c_short)
    !> The size in bytes of an integer 
    integer(c_size_t), parameter :: bytes_int             = c_sizeof(1_c_int)
    !> The size in bytes of a long integer
    integer(c_size_t), parameter :: bytes_long_int        = c_sizeof(1_c_long)
    !> The size in bytes of single precision real
    integer(c_size_t), parameter :: bytes_single_real     =  c_sizeof(1.0_c_float)
    !> The size in bytes of double precision real
    integer(c_size_t), parameter :: bytes_double_real     =  c_sizeof(1.0_c_double)
    !> The size in bytes of single precision complex
    integer(c_size_t), parameter :: bytes_single_complex  =  2_c_size_t * bytes_single_real
    !> The size in bytes of double precision complex
    integer(c_size_t), parameter :: bytes_double_complex  =  2_c_size_t * bytes_double_real


interface 

    !> Allocates a memory and returns C pointer to that memory
    !> @param[out] memory - the allocated memory
    !> @param[in]  memsize - the size in bytes to allocate
    !> @param[in]  device_id - the device id in which the allocation will happen
    module subroutine allocate_device_memory(memory, memsize, device_id)
        type(c_ptr), intent(out) :: memory 
        integer(c_size_t), intent(in) :: memsize
        integer, intent(in) :: device_id
    end subroutine allocate_device_memory

    !> Deallocates a memory pointed by 
    !> @param[inout] memory - the allocated memory
    !> @param[in]  device_id - the device id in which the allocation will happen
    module subroutine deallocate_device_memory(memory, device_id)
        type(c_ptr), intent(inout) :: memory 
        integer, intent(in) :: device_id
    end subroutine deallocate_device_memory

    !> Given data in the host it provides the C_ptr in the device for device backend
    !> otherwise it is equivalent to c_loc
    !> @param[in]  host_data - the allocated memory in host for which the associated device ptr is searched
    !> @param[in]  device_id - the device id in which the allocation will happen
    !> @result device_c_ptr to the mapped data in the device or c_loc(host_data) in the CPU backend
    module type(c_ptr) function get_device_pointer(host_data, device_id)
        type(*), target, intent(in) :: host_data(..)
        integer, intent(in) :: device_id
    end function get_device_pointer


    !  @brief Fortran interface for C++ function `generate_batched_array`.
    !
    !  This subroutine serves as a Fortran wrapper for the C++ function
    !  `generate_batched_array`, which generates pointers to batched data
    !  segments stored in a contiguous memory block. The batched data layout
    !  is assumed to have the batch identifier as the slowest varying index.
    !
    !  @param[in] data            Base address of the contiguous data block.
    !  @param[in] nbatch          Number of batches.
    !  @param[in] batch_byte_size Size (in bytes) of each batch.
    !  @param[out] batched_ptr    Array of pointers where each entry points
    !                             to the start of a batch within `data`.
    ! -----------------------------------------------------------------------------
    subroutine generate_batched_array(data, nbatch, batch_byte_size, batched_ptr) &
        bind(C, name="generate_batched_array")
        use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
        implicit none
        ! Arguments
        type(c_ptr), value       :: data              ! Base address of the data
        integer(c_size_t), value :: nbatch            ! Number of batches
        integer(c_size_t), value :: batch_byte_size   ! Byte size of each batch
        type(c_ptr)              :: batched_ptr(*)    ! Array of batch pointers
    end subroutine generate_batched_array

end interface

! Interfaces to C routines
interface
    ! Interface for malloc function
    function c_malloc(size) bind(C, name="malloc")
        import :: c_size_t, c_ptr
        type(c_ptr) :: c_malloc
        integer(c_size_t), value :: size
    end function c_malloc

    ! Interface for free function
    subroutine c_free(ptr) bind(C, name="free")
        import :: c_ptr
        type(c_ptr), value :: ptr
    end subroutine c_free
end interface

contains

    !> Returns the size in bytes for an arbitrary kind
    !> @param[in] data - the data from which the size in bytes is retrieved 
    !> @return : bytes in memory of data object
    IFORT_BUG_SAFE_PURE integer(c_size_t) function bytes_size(data)
        IFORT_BUG_POLYCLASS_ASSUMED_RANK, intent(in) :: data(..)

        select rank(data)
            rank(0)
                bytes_size = storage_size(data, kind=c_size_t) / bit_size
            rank default
                bytes_size = size(data) * storage_size(data,kind=c_size_t) / bit_size
        end select

    end function bytes_size


    function get_pointer_alignment(c_pointer) result(alignment)
        type(c_ptr), intent(in) :: c_pointer

        integer(c_int) :: alignment

        integer(c_intptr_t) :: pointer_int

        ! Transform the c_ptr to an integer
        pointer_int = transfer(c_pointer, pointer_int)

        if      (modulo(pointer_int,64) == 0) then
            alignment = 64
        else if (modulo(pointer_int,32) == 0) then
            alignment = 32
        else if (modulo(pointer_int,16) == 0) then
            alignment = 16
        else if (modulo(pointer_int,8) == 0)  then
            alignment = 8
        else if (modulo(pointer_int,4) == 0)  then
            alignment = 4
        else if (modulo(pointer_int,2) == 0)  then
            alignment = 2
        else
            error stop "get_pointer_alignment : the pointer is not aligned"
        end if
        return

    end function get_pointer_alignment

end module m_memory_device
