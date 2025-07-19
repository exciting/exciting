! Copyright (C) 2020-2024 GreenX library
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
!> This file contains a module for exposing the idiel_t to C
!> the aim is to create a C++ library and a Python module using the C interface

!> This module defines the subroutines for exposing the idiel_t to C
!> Notice that names should be lower case as the compiled functions will be named using lower case
!> as of Fortran ISO interfacing with C
module idiel_f90_C

    use iso_c_binding
    use idiel_constants, only: i32, r64
    use idiel, only: idiel_t
    
    implicit none
    
    public
    
contains
    
    !> Allocates idiel_t
    subroutine allocate_idiel_t(object_ptr) bind(C)
        type(c_ptr), intent(out) :: object_ptr
        type(idiel_t), pointer :: object
        
        allocate(object)
        object_ptr = C_loc(object)

    end subroutine allocate_idiel_t
    
    !> Deallocates the idiel_t
    subroutine deallocate_idiel_t(object_ptr) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        type(idiel_t), pointer     :: object
        
        call C_F_pointer(object_ptr, object)
        call object%clean()

    end subroutine deallocate_idiel_t
    
    subroutine init_common(object_ptr, lattice, natoms, redpos, elements, nq, dim) bind(C)
        type(c_ptr),  intent(inout)      :: object_ptr
        real(r64),    intent(in)         :: lattice(3,3)
        integer(i32), intent(in), value  :: natoms
        real(r64),    intent(in)         :: redpos(3,natoms) 
        integer(i32), intent(in)         :: elements(natoms)
        integer(i32), intent(in)         :: nq(3)
        integer(i32), intent(in), value  :: dim
        
        type(idiel_t), pointer :: object

        call C_F_pointer(object_ptr, object)
                
        call object%init_common(lattice, redpos, elements, nq, dim)

    end subroutine init_common
    
    subroutine set_dielectric_blocks(object_ptr, h, wl, wu, ib) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        complex(r64), target, intent(in)                :: h(:,:)
        complex(r64), target, intent(in)                :: wl(:,:)
        complex(r64), target, intent(in)                :: wu(:,:)
        complex(r64), target, optional, intent(in)      :: ib(:,:)
        
        type(idiel_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        if (present(ib)) then
            call object%set_dielectric_blocks(h, wl, wu, ib)
        else 
            call object%set_dielectric_blocks(h, wl, wu)
        end if 

    end subroutine set_dielectric_blocks

    subroutine compute_anisotropic_avg_inversedielectric_3d(object_ptr, hermitian) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        logical,     intent(in)    :: hermitian
        
        type(idiel_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%compute_anisotropic_avg_inversedielectric_3d(hermitian)
        
    end subroutine compute_anisotropic_avg_inversedielectric_3d

    subroutine compute_anisotropic_avg_scrcoulomb_2d(object_ptr, hermitian) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        logical,     intent(in)    :: hermitian
        
        type(idiel_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%compute_anisotropic_avg_scrcoulomb_2d(hermitian)
        
    end subroutine compute_anisotropic_avg_scrcoulomb_2d
    
    subroutine invert_body(object_ptr, body) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        complex(r64), allocatable, intent(in) :: body(:,:)
        
        type(idiel_t), pointer :: object
        call C_F_pointer(object_ptr, object)

        call object%invert_body(body)

    end subroutine invert_body

    function get_n_basis(object_ptr) result(nbasis) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        integer(i32) :: nbasis
        
        type(idiel_t), pointer :: object
        call C_F_pointer(object_ptr, object)

        nbasis = object%get_n_basis()

    end function get_n_basis

    function head(object_ptr) result(ptr) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        type(c_ptr) :: ptr

        type(idiel_t), pointer :: object
        
        call C_F_pointer(object_ptr, object)
        ptr = C_loc(object%idiel_head)

    end function head

    function wing_lower(object_ptr) result(ptr) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        type(c_ptr) :: ptr

        type(idiel_t), pointer :: object

        call C_F_pointer(object_ptr, object)
        ptr = C_loc(object%idiel_wingL)

    end function wing_lower

    function wing_upper(object_ptr) result(ptr) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        type(c_ptr) :: ptr

        type(idiel_t), pointer :: object

        call C_F_pointer(object_ptr, object)
        ptr = C_loc(object%idiel_wingU)

    end function wing_upper

    function body(object_ptr) result(ptr) bind(C)
        type(c_ptr), intent(inout) :: object_ptr
        type(c_ptr) :: ptr

        type(idiel_t), pointer :: object

        call C_F_pointer(object_ptr, object)
        ptr = C_loc(object%idiel_body)

    end function body

end module idiel_f90_C
