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
!> This module contains information regarding crystal symmetry agnostic of any DFT solver
module idiel_crystal_symmetry    

    use idiel_constants,    only: i32, aip, twopi, zzero
    use iso_c_binding,  only: C_int, C_double, C_char
    use idiel_crystal_cell, only: cell_t

    implicit none

    private
    public  symmetry_t

    !> This class contains the symmetry needed by the programm
    type symmetry_t
        !> Space group id
        integer(i32) :: spgid
        !> Number of symmetry operations
        integer(i32) :: nsym
        !> Rotation matrices crystal coordinates
        integer(i32), allocatable :: rot(:,:,:)
        !> Rotation matrices for qpoints (no inversion symmetry)
        real(aip), allocatable :: qrot(:,:,:)
        !> Rotation matrices Cartesian coordinates
        real(aip), allocatable :: crot(:,:,:)
    contains
#ifdef USE_SPGLIB
        procedure, public :: initialize=>init_symmetry
#endif
        procedure, public :: symmetryze_complex_tensor
        final :: clean
    end type symmetry_t

#ifdef USE_SPGLIB

#ifdef USE_SINGLE_PRECISION
#define _gesv sgesv
#else
#define _gesv dgesv
#endif

interface

    !> Interface to spg_get_multiplicity (this gives the number of symmetry operations)
    !> @param[in] lattice - crystal lattice
    !> @param[in] positions - reduced positions
    !> @param[in] types     - atomic types
    !> @param[in] natoms    - number of atoms
    !> @param[in] symprec   - tolerance for the symmetry search
    !> @returns   symmetry_get_multiplicity - number of symmetry operations
    function symmetry_get_multiplicity(lattice, positions, types, natoms, symprec) bind (C, &
                                       name="spg_get_multiplicity")
        import 
        real(C_double)        :: lattice(3,3)
        real(C_double)        :: positions(3,natoms)
        integer(C_int)        :: types(natoms)
        integer(C_int), value :: natoms
        real(C_double), value :: symprec
        integer(C_int)        :: symmetry_get_multiplicity
    end function symmetry_get_multiplicity

    !> Interface to spg_get_symmetry which gives the symmetry operations
    !> @param[in,out] rotations  - rotations in crystal coordinates
    !> @param[in,out] translation - fractional translations in crystal coordinates
    !> @param[in] nops      - number of symmetry operations
    !> @param[in] lattice - crystal lattice
    !> @param[in] positions - reduced positions
    !> @param[in] types     - atomic types
    !> @param[in] natoms    - number of atoms
    !> @param[in] symprec   - tolerance for the symmetry search
    !> @returns   symmetry_compute - number of symmetry operations
    function symmetry_compute(rotations, translations, nops, lattice, positions, types, natoms,  & 
                              symprec) bind (C,name="spg_get_symmetry")
        import 
        integer(C_int)        :: rotations(3,3,nops)
        real(C_double)        :: translations(3,nops)
        integer(C_int), value :: nops
        real(C_double)        :: lattice(3,3)
        real(C_double)        :: positions(3,natoms)
        integer(C_int)        :: types(natoms)
        integer(C_int), value :: natoms
        real(C_double), value :: symprec
        integer(C_int)        :: symmetry_compute
    end function symmetry_compute

    !> Interface to spg_get_international which gives the symmetry operations
    !> @param[out] symbol   - the spacegroup symbol
    !> @param[in] lattice   - crystal lattice
    !> @param[in] positions - reduced positions
    !> @param[in] types     - atomic types
    !> @param[in] natoms    - number of atoms
    !> @param[in] symprec   - tolerance for the symmetry search
    !> @returns   get_spacegroup - space group number
    function get_spacegroup(symbol, lattice, positions, types, natoms, &
                            symprec) bind (C, name="spg_get_international")
        import
        integer(C_int), value  :: natoms
        character(C_char)      :: symbol(11)
        real(C_double)         :: lattice(3,3)
        real(C_double)         :: positions(3,natoms)
        integer(C_int)         :: types(natoms)
        real(C_double), value  :: symprec
        integer(C_int)         :: get_spacegroup
   end function get_spacegroup


end interface
#endif

contains

    !> @brief Cleans the elements in the cell type
    !> @param[in,out] this - the symmetry object to deallocate
    subroutine clean(this)
        type(symmetry_t) :: this
        if (allocated(this%rot)) deallocate(this%rot)
        if (allocated(this%crot)) deallocate(this%crot)
        if (allocated(this%qrot)) deallocate(this%qrot)
    end subroutine clean

#ifdef USE_SPGLIB

    !> @brief Initializes the symmetry object
    !> @param[in,out] this      - the cell object to deallocate
    !> @param[in] cell          - cell object from which the symmetry should be searched
    !> @param[in] tolerance - the tolerance to determine the symmetry of the cell (optional)
    subroutine init_symmetry(this, cell, tolerance)
        
        class(symmetry_t), intent(inout) :: this
        type(cell_t), intent(in) :: cell
        real(aip), intent(in), optional :: tolerance

        !Locals
        real(aip) :: tol = 1.0e-5_aip
        integer(i32) :: isym
        integer(C_int), allocatable ::  rot(:,:,:)
        real(C_double), allocatable ::  tau(:,:)
        character(C_char) :: symbol(11) 
        real(aip) :: lat(3,3), invlat(3,3)
        real(aip) :: rT(3,3), temp(3,3), hermitian(3,3)
        integer :: info, ipiv(3)

        ! External
        external :: _gesv

        ! If present change the tolerance to the desired value
        if (present(tolerance))  tol = tolerance

        this%spgid = int(get_spacegroup(symbol,                        &
                                        real(cell%lattice,C_double) ,  &
                                        real(cell%redpos,C_double)  ,  &
                                        int(cell%elements,C_int),      &   
                                        int(cell%natoms,C_int),        & 
                                        real(tol,C_double)), i32)



        this%nsym = int(symmetry_get_multiplicity(real(cell%lattice,C_double) ,  &
                                                  real(cell%redpos,C_double)  ,  &
                                                  int(cell%elements,C_int),      &   
                                                  int(cell%natoms,C_int),        & 
                                                  real(tol,C_double)), i32)
        
        allocate(rot(3,3,this%nsym), tau(3,this%nsym))

        this%nsym = int(symmetry_compute(rot,                           &
                                         tau,                           &
                                         int(this%nsym,C_int),          &
                                         real(cell%lattice,C_double) ,  &
                                         real(cell%redpos,C_double)  ,  &
                                         int(cell%elements,C_int),      &   
                                         int(cell%natoms,C_int),        & 
                                         real(tol,C_double)), i32)
        
        ! Give it in Fortran Order
        do isym = 1, this%nsym
            rot(:,:,isym) = transpose(rot(:,:,isym))
        end do

        allocate(this%rot(3,3,this%nsym))
        allocate(this%qrot(3,3,this%nsym))
        this%rot = int(rot,i32)

        ! Compute the rotation in the Cartesian rotations
        allocate(this%crot(3,3,this%nsym))
        
        lat    = transpose(cell%lattice)
        invlat = transpose(cell%rlattice / twopi)
        
        do isym = 1, this%nsym
            this%crot(:,:,isym) = transpose(matmul(lat,this%rot(:,:,isym)))
            this%crot(:,:,isym) = transpose(matmul(invlat,this%crot(:,:,isym)))
        end do

        ! Compute the qrotations (crystal coordinates)
        hermitian = matmul(cell%lattice,transpose(cell%lattice))
        do isym = 1, this%nsym
            rT = transpose(real(this%rot(:,:,isym),aip))
            temp = hermitian
            call _gesv(3,3,temp,3,ipiv,rT,3,info)
            this%qrot(:,:,isym) = transpose(matmul(rT,hermitian))
        end do

    end subroutine init_symmetry

#endif

    !> Symmetrizes a Cartesian complex double precision tensor
    !> @param[in] this - symmetry object
    !> @param[in] mat33 - the matrix to symmetrize
    !> @result syidiel_mat33 - the symmetrized vector
    function symmetryze_complex_tensor(this, mat33) result(sym_mat33)
        class(symmetry_t), intent(in) :: this
        complex(aip), intent(in) :: mat33(3,3)
        
        complex(aip) :: sym_mat33(3,3)
        
        integer(i32) :: isym

        sym_mat33 = zzero

        do isym = 1, this%nsym
            sym_mat33 = sym_mat33 + matmul(transpose(this%crot(:,:,isym)),matmul(mat33,this%crot(:,:,isym))) 
        end do

        sym_mat33 = sym_mat33 / this%nsym

    end function symmetryze_complex_tensor


end module idiel_crystal_symmetry
