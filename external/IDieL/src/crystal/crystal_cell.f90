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
!> This module contains modules for crystalline structure

!> Module containing the information about crystal structure except for symmetry
module idiel_crystal_cell

    use idiel_constants, only: i32, aip, twopi

    implicit none

    private 
    public  cell_t

    !> This class contains the structural information 
    type cell_t 
        !> The lattice each row is a lattice vector (Fortran order)
        real(aip) :: lattice(3,3)
        !>  The reciprocal lattice (Fortran order)
        real(aip) :: rlattice(3,3) 
        !>  The lattice volume
        real(aip) :: vuc
        !>  Area between b1xb2 vectors (reciprocal)
        real(aip) :: area_r_12
        !>  Area between b1xb3 vectors (reciprocal)
        real(aip) :: area_r_13
        !>  Area between b2xb3 vectors (reciprocal)
        real(aip) :: area_r_23
        !>  The number of atoms in the unit cell
        integer(i32) :: natoms
        !>  The reduced positions (3,natoms),  this in Fortran is already compatible with what expected from spglib
        real(aip), allocatable :: redpos(:,:) 
        !> The id of elements
        integer(i32), allocatable :: elements(:)
    contains
        procedure, public :: initialize=>init_cell_t, get_kmax_subcell_bz
        final :: clean
    end type cell_t


contains

    !> Cleans the elements in the cell_t type
    !> @param[in,out] this - the cell_t object to deallocate
    subroutine clean(this)

        type(cell_t), intent(inout) :: this

        if(allocated(this%redpos)) deallocate(this%redpos)
        if(allocated(this%elements)) deallocate(this%elements)

    end subroutine clean

    !> Initializes the cell_t object
    !> @param[in,out] this      - the cell_t object to deallocate
    !> @param[in]     lattice_  - lattice vectors given in rows [nm]
    !> @param[in]     redpos_   - reduced positions (3,natoms)
    !> @param[in]     elements_ - list of elements
    subroutine init_cell_t(this, lattice_, redpos_, elements_)
        class(cell_t), intent(inout) :: this
        real(aip), intent(in)        :: lattice_(3,3)
        real(aip), intent(in)        :: redpos_(:,:) 
        integer(i32), intent(in)     :: elements_(:)

        ! Locals
        real(aip) :: a(3), b(3), c(3)
        real(aip) :: ap(3), bp(3), cp(3)

        ! Init class variables
        this%lattice  = lattice_
        allocate(this%redpos, source = redpos_)
        allocate(this%elements, source = elements_)
        this%natoms   = size(this%redpos,2)

        ! The elements of the lattice
        a(:) = lattice_(1,:)
        b(:) = lattice_(2,:)
        c(:) = lattice_(3,:)

        ! Compute the volume
        this%vuc = dot_product(a,cross(b,c))

        !Compute the reciprocal lattice, volumes, and related information
        ap = twopi * cross(b,c) / this%vuc
        bp = twopi * cross(c,a) / this%vuc
        cp = twopi * cross(a,b) / this%vuc

        this%area_r_12 = norm2(cross(ap,bp))
        this%area_r_13 = norm2(cross(ap,cp))
        this%area_r_23 = norm2(cross(bp,cp))

        this%vuc = abs(this%vuc)

        this%rlattice(1,:) = ap(:)
        this%rlattice(2,:) = bp(:)
        this%rlattice(3,:) = cp(:)

    contains 
        
        !> Cross product function
        !> @param[in] v1 - the first 3d vector
        !> @param[in] v2 - the second 3d vector
        !> @result    v3 - the v1 X v2
        function cross(v1,v2) result(v3)
            
            real(aip), intent(in) :: v1(3)
            real(aip), intent(in) :: v2(3)
            real(aip) :: v3(3)

            v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
            v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
            v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
            
        end function cross 


    end subroutine init_cell_t


    !> Given a subcell in the first BZ get the distance from the cell of the subcell to the surface
    !> @param[in]     this      - the cell_t object to deallocate
    !> @param[in]     nqgrid    - the size of the regular BZ-mesh
    !> @param[in]     kdir      - the direction in which to look for the boundary [Cartesian]
    !> @param[out]    kmax      - the distance to the subcell surface     
    pure subroutine get_kmax_subcell_bz(this, nqgrid, kdir, kmax)
        
        class(cell_t), intent(in)           :: this
        integer(i32),  intent(in)           :: nqgrid(3)
        real(aip),     intent(in)           :: kdir(:,:)
        real(aip), allocatable, intent(out) :: kmax(:)

        ! Locals
        real(aip), allocatable :: kred(:,:) ! k in reduced coordinates (3,nr)
        integer(i32) :: nr                  ! number of points
        integer(i32) :: i, j
        real(aip) :: lambda_candidates(3), lambda ! The scaling factor in the line equation, i.e. (kx,ky,kz) = (0,0,0) + lambda * kdir  
        real(aip) :: k_check              ! The value to check for the intersection with one of the planes 

        ! Get the number of points
        nr = size(kdir,1) 

        ! Allocate space for kmax
        allocate(kmax(nr))

        ! First compute the xyz points in reciprocal units
        kred = matmul(this%lattice / twopi , transpose(kdir) )
        
        ! Iterate over each point
        do i = 1, nr
            ! Iterate over the reciprocal axis
            do j = 1, 3
                ! Compute all lambdas that provides a collision with one of the planes defining the boundaries
                ! Handle vectors parallel to the planes
                if (abs(kred(j,i)) .lt. 1.0e-12) then 
                    kred(j,i) = 0.0_aip
                    lambda_candidates(j) = huge(1.0_aip)
                else
                    ! The collision points for the planes defining a subcell (given uniform division) of the unit cube centered in 0 is by definition sgn(kdir(i)) * 0.5/nq(i)
                    ! The sign is to constrain our search to lambda > 0
                    k_check = sign(0.5_aip / nqgrid(j), kred(j,i)) 
                    lambda_candidates(j) = k_check / kred(j,i)
                end if
            end do 
            ! The minum lamba will provide the point in the subcell surface 
            lambda = minval(lambda_candidates)
            ! Return to reduced point to Cartesian units and compute the distance 
            kmax(i) = norm2(lambda * kdir(i,:))
        end do

        ! The top algorithm can be rewritten  as:
        ! do i = 1, nr
        !     kmax(i) = 1.0/ (2.0 * maxval(abs(nqgrid(:) * kred(:,i)), 1))
        ! end do 

    end subroutine get_kmax_subcell_bz

end module idiel_crystal_cell
