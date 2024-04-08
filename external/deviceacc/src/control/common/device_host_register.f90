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
!> Contains Fortran interface for C functions contained in device_host_register.cpp
!> See device_host_register.cpp for the detailed documentation
module m_device_host_register_fortran

    use iso_c_binding

    implicit none

    private
    public  :: device_host_register

    type device_host_register
        !> The C ptr holding the pointer to cxx class
        type(c_ptr) :: device_host_register_cxx = c_null_ptr
    contains
        procedure, public :: init, finish
        procedure, public :: alloc, remove
        procedure, public :: assoc, deassoc
        procedure, public :: to_device, from_device
        procedure, public :: device_ptr
    end type device_host_register

interface
    subroutine constructor_device_host_register(rg) bind(c, name="_constructor_device_host_register")
        import
        type(c_ptr) :: rg
    end subroutine constructor_device_host_register

    subroutine destructor_device_host_register(rg) bind(c, name="_destructor_device_host_register")
        import
        type(c_ptr), value :: rg
    end subroutine destructor_device_host_register

    subroutine alloc_device_device_host_register(rg, id_len, id, size, device_id) bind(c, name="_alloc_device_device_host_register")
        import
        type(c_ptr), value                               :: rg
        integer(c_int) , value                           :: id_len
        character(kind=c_char), dimension(*), intent(in) :: id
        integer(c_size_t), value                         :: size
        integer(c_int), value                            :: device_id
    end subroutine alloc_device_device_host_register

    subroutine associate_device_device_host_register(rg, id_len, id, host_ptr) bind(c, name="_associate_device_device_host_register")
        import
        type(c_ptr), value                               :: rg
        integer(c_int) , value                           :: id_len
        character(kind=c_char), dimension(*), intent(in) :: id
        type(c_ptr), value                               :: host_ptr
    end subroutine associate_device_device_host_register

    subroutine disassociate_device_device_host_register(rg, id_len, id) bind(c, name="_disassociate_device_device_host_register")
        import
        type(c_ptr), value                               :: rg
        integer(c_int) , value                           :: id_len
        character(kind=c_char), dimension(*), intent(in) :: id
    end subroutine disassociate_device_device_host_register

    subroutine host_to_device_device_device_host_register(rg, id_len, id) bind(c, name="_host_to_device_device_device_host_register")
        import
        type(c_ptr), value                               :: rg
        integer(c_int) , value                           :: id_len
        character(kind=c_char), dimension(*), intent(in) :: id
    end subroutine host_to_device_device_device_host_register

    subroutine device_to_host_device_device_host_register(rg, id_len, id) bind(c, name="_device_to_host_device_device_host_register")
        import
        type(c_ptr), value                               :: rg
        integer(c_int) , value                           :: id_len
        character(kind=c_char), dimension(*), intent(in) :: id
    end subroutine device_to_host_device_device_host_register

    subroutine remove_device_device_host_register(rg, id_len, id) bind(c, name="_remove_device_device_host_register")
        import
        type(c_ptr), value                               :: rg
        integer(c_int) , value                           :: id_len
        character(kind=c_char), dimension(*), intent(in) :: id
    end subroutine remove_device_device_host_register

    function get_device_ptr_device_device_host_register(rg, id_len, id) bind(c, name="_get_device_ptr_device_device_host_register")
        import
        type(c_ptr), value                               :: rg
        integer(c_int) , value                           :: id_len
        character(kind=c_char), dimension(*), intent(in) :: id
        type(c_ptr) :: get_device_ptr_device_device_host_register
    end function get_device_ptr_device_device_host_register

end interface

contains

    subroutine init(this)
        class(device_host_register), intent(inout) :: this
        call constructor_device_host_register(this%device_host_register_cxx)
    end subroutine init

    subroutine finish(this)
        class(device_host_register), intent(inout) :: this
        call destructor_device_host_register(this%device_host_register_cxx)
        this%device_host_register_cxx = c_null_ptr
    end subroutine finish

    subroutine alloc(this, id, size, device_id)
        class(device_host_register), intent(inout) :: this
        character(len=*), intent(in) :: id
        integer(c_size_t), value :: size
        integer(c_int), value :: device_id
        call alloc_device_device_host_register(this%device_host_register_cxx, len(id, kind=c_int), id, size, device_id)
    end subroutine alloc

    subroutine remove(this, id)
        class(device_host_register), intent(inout) :: this
        character(len=*), intent(in) :: id
        call remove_device_device_host_register(this%device_host_register_cxx, len(id, kind=c_int), id)
    end subroutine remove

    subroutine assoc(this, id, host_ptr)
        class(device_host_register), intent(inout) :: this
        character(len=*), intent(in) :: id
        type(c_ptr), value :: host_ptr
        call associate_device_device_host_register(this%device_host_register_cxx, len(id, kind=c_int), id, host_ptr)
    end subroutine assoc

    subroutine deassoc(this, id)
        class(device_host_register), intent(inout) :: this
        character(len=*), intent(in) :: id
        call disassociate_device_device_host_register(this%device_host_register_cxx, len(id, kind=c_int), id)
    end subroutine deassoc

    subroutine to_device(this, id)
        class(device_host_register), intent(inout) :: this
        character(len=*), intent(in) :: id
        call host_to_device_device_device_host_register(this%device_host_register_cxx, len(id, kind=c_int), id)
    end subroutine to_device

    subroutine from_device(this, id)
        class(device_host_register), intent(inout) :: this
        character(len=*), intent(in) :: id
        call device_to_host_device_device_host_register(this%device_host_register_cxx, len(id, kind=c_int), id)
    end subroutine from_device

    function device_ptr(this, id) result(d_ptr)
        class(device_host_register), intent(inout) :: this
        character(len=*), intent(in) :: id
        type(c_ptr) :: d_ptr
        d_ptr = get_device_ptr_device_device_host_register(this%device_host_register_cxx, len(id, kind=c_int), id)
    end function device_ptr

end module m_device_host_register_fortran
