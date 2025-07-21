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
!> Contains procedures for averages of the dielectric matrix

!> This module contains the procedures for the computation of the dielectric matrix averages
module idiel

    use idiel_constants, only: i32, r64, aip, pi, onethird, twopi, fourpi, zzero, zone, iunit
    use idiel_crystal_cell, only: cell_t
    use idiel_crystal_symmetry, only: symmetry_t
    use idiel_sph_quadrature, only: compute_angular_mesh_lebedev_131, compute_angular_mesh_lebedev_41
    use idiel_spherical_harmonics, only: sph_harm, sph_harm_expansion
    use idiel_circle_quadrature, only: compute_angular_mesh_gauss_legendre
    use idiel_circular_harmonics, only: circ_harm, circ_harm_expansion
    use m_device_world_t, only: device_world_t

    implicit none

    private
    public idiel_t


    !> Order of harmonic spherical/circular (Fourier) expansion
    integer(i32), parameter :: lmax = 20_i32
    !> Number of spherical harmonics
    integer(i32), parameter :: nsph = (lmax + 1)**2
    !> Number of pairs spherical harmonics
    integer(i32), parameter :: nsph_pair = 231_i32
    !> Number of circular harmonics
    integer(i32), parameter :: ncir = (1 + 2*lmax)
    !> Number of radial points for 2D integrals
    integer(i32), parameter :: nr   = 151_i32 
    !> Number of points for small 2D circular (i.e. Fourier) expansion
    integer(i32), parameter :: size_mesh_2d_coarse = 75_i32
    !> Number of points for the big 2D circular (i.e. Fourier) expansion
    integer(i32), parameter :: size_mesh_2d_fine   = 255_i32

    !> This type is used to perform all the IDieL operations
    type idiel_t
        !> Cell information
        type(cell_t), private :: cell
        !> Symmetry information
        type(symmetry_t), private :: symmetry
        !> The BZ mesh size
        integer(i32), private :: nq(3) 
        !> Angular mesh
        real(aip), allocatable, private :: ang(:,:) 
        !> Mesh points (Cartesian)
        complex(aip), allocatable, private :: xyz(:,:) 
        !> Weights for the integrals (fine mesh)
        real(aip), allocatable, private :: weights_fine(:)
        !> Weights for the integrals (coarse mesh)
        real(aip), allocatable, private :: weights(:)
        !> Spherical/circular harmonics
        complex(aip), allocatable, private :: ylm(:,:)
        !> 3D angular integrals
        complex(aip), allocatable, private :: angular_integrals(:)
        !> Pointer to the head of the dielectric matrix
        complex(aip), pointer, private :: head(:,:)    => null()
        !> Pointer to the lower wing, i.e. in PW G0
        complex(aip), pointer, private :: wingL(:,:)   => null()
        !> Pointer to the upper wing, i.e. in PW 0G 
        complex(aip), pointer, private :: wingU(:,:)   => null()
        !> Pointer to the inverse of the body block of the dielectric matrix 
        complex(aip), pointer, private :: Binv(:,:)    => null()
        !> Actual data in case that it has been computed using this library
        complex(aip), allocatable :: Binv_data(:,:)
        !> The head of the inverse
        complex(aip) :: idiel_head
        !> The lower wing of the inverse dielectric matrix / screened Coulomb
        complex(aip), allocatable :: idiel_wingL(:)
        !> The upper wing of the inverse dielectric matrix / screened Coulomb
        complex(aip), allocatable :: idiel_wingU(:)
        !> The body of the inverse dielectric matrix / screened Coulomb
        complex(aip), allocatable :: idiel_body(:,:)
        !> Number of points of the quadrature
        integer(i32), private :: quadrature_npoints 
        !> Pointer to the handler of GPU-CPU
        type(device_world_t), pointer, private :: world
        !> Logical to know if the world is external or allocated in the object
        logical, private :: internal_world
        !> Dimensionaly
        integer(i32), private :: dim = 3
        !> 2D angle for large integral
        real(aip), private, allocatable :: phi(:)
        !> 2D radial mesh for interpolation
        real(aip), private, allocatable :: radii(:)
        !> 2D rmax as function of angle
        real(aip), allocatable, private :: rmax2d(:)
        !> rcut for the 2D Coulomb potential
        real(aip), private :: rcut 
        !> Exponential function of the screened Coulomb potential 
        complex(aip), allocatable, private :: vr(:,:)
        !> Circular harmonics in the coarse mesh
        complex(aip), allocatable, private :: blm_coarse(:,:)
        !> Circular harmonics in the fine mesh
        complex(aip), allocatable, private :: blm_fine(:,:)
    contains
        procedure, public  :: init_common, set_dielectric_blocks, invert_body, get_n_basis
        procedure, public  :: compute_anisotropic_avg_inversedielectric_3d
        procedure, private :: compute_anisotropic_avg_hermitian, compute_anisotropic_avg_general
        procedure, public  :: compute_anisotropic_avg_scrcoulomb_2d
        procedure, private :: compute_anisotropic_avg_scrcoulomb_2d_general, compute_anisotropic_avg_scrcoulomb_2d_hermitian
        procedure, public  :: clean
    end type idiel_t

interface

    !> This initializes everything common (i.e. without frequency dependency) to all computations
    !> @param[in,out] this - idiel_t object
    !> @param[in]     lattice       - lattice vectors given in rows
    !> @param[in]     redpos        - reduced positions (3,natoms)
    !> @param[in]     elements      - list of elements
    !> @param[in]     nq            - the BZ mesh
    !> @param[in]     dim           - the dim of the system
    !> @param[in]     nysm          - the number of symmetry operations
    !> @param[in]     crot          - the rotations in Cartesian coordinates (last index is the space group operation id)
    !> @param[in]     device_world  - the GPU-CPU communicator
    !> @param[in]     host_world    - the MPI communicator
    module subroutine init_common(this, lattice, redpos, elements, nq, dim, nsym, crot, device_world, host_world)
        use iso_c_binding, only: c_int
        class(idiel_t), target, intent(inout)              :: this
        real(aip), intent(in)                              :: lattice(3,3)
        real(aip),  intent(in)                             :: redpos(:,:)
        integer(i32),  intent(in)                          :: elements(:)
        integer(i32), intent(in)                           :: nq(3)
        integer(i32), intent(in), optional                 :: dim
        integer(i32), intent(in), optional                 :: nsym
        real(aip), intent(in), optional                    :: crot(:,:,:)
        type(device_world_t), intent(in), target, optional :: device_world
        integer(c_int), intent(in), optional               :: host_world
    end subroutine init_common

    !> This nullify and deallocates the objects
    !> @param[in] this - idiel_t object
    module subroutine clean(this)
        class(idiel_t), intent(inout) :: this
    end subroutine clean

    !> This subroutine inits pointers to a given value
    !> @param[in] this - idiel_t object
    !> @param[in] h    - head of the dielectric matrix
    !> @param[in] wl   - lower wing of the dielectric matrix
    !> @param[in] wu   - upper wing of the dielectric matrix
    !> @param[in] ib   - inverse of the body of the dielectric matrix (optional)
    module subroutine set_dielectric_blocks(this, h, wl, wu, ib)
        class(idiel_t), intent(inout)              :: this
        complex(aip), target, intent(in)           :: h(:,:)
        complex(aip), target, intent(in)           :: wl(:,:)
        complex(aip), target, intent(in)           :: wu(:,:)
        complex(aip), target, optional, intent(in) :: ib(:,:)
    end subroutine set_dielectric_blocks

    !> This computes the average of the anisotropic inverse dielectric matrix around Gamma
    !> @param[in] this       - the current idiel_t object for which to compute the average
    !> @param[in] hermitian  - is the dielectric matrix hermitian
    module subroutine compute_anisotropic_avg_inversedielectric_3d(this, hermitian)
        class(idiel_t), intent(inout) :: this
        logical,           intent(in) :: hermitian
    end subroutine compute_anisotropic_avg_inversedielectric_3d
    
    !> This computes the average of the anisotropic inverse dielectric matrix around Gamma for Hermitian matrices (i.e. purely imaginary frequencies)
    !> @param[in] this       - the current idiel_t object for which to compute the average
    module subroutine compute_anisotropic_avg_hermitian(this)
        class(idiel_t), target, intent(inout) :: this
    end subroutine compute_anisotropic_avg_hermitian


    !> This computes the average of the anisotropic inverse dielectric matrix around Gamma for general matrices 
    !> @param[in] this       - the current idiel_t object for which to compute the average
    module subroutine compute_anisotropic_avg_general(this)
        class(idiel_t), target, intent(inout) :: this
    end subroutine compute_anisotropic_avg_general

    !> Inverts the body and stores it (GPU or CPU depending on the compilation)
    !> @param[in] this - the idiel_t in which to store the inverse of body
    !> @param[in] body - the body to invert
    module subroutine invert_body(this, body)
        class(idiel_t), intent(inout), target :: this
        complex(aip), intent(in) :: body(:,:)
    end subroutine invert_body

    !> This function returns the basis size (for CXX and Python compatibility)
    !> @param[in] this - the idiel_t object to check the size of the used basis
    module function get_n_basis(this) result(nbasis)
        class(idiel_t), intent(inout), target :: this
        integer(i32) :: nbasis
    end function

    module subroutine compute_anisotropic_avg_scrcoulomb_2d(this, hermitian)
        class(idiel_t), intent(inout) :: this
        logical, intent(in)           :: hermitian
    end subroutine compute_anisotropic_avg_scrcoulomb_2d

    module subroutine compute_anisotropic_avg_scrcoulomb_2d_general(this)
        class(idiel_t), target, intent(inout) :: this
    end subroutine compute_anisotropic_avg_scrcoulomb_2d_general

    module subroutine compute_anisotropic_avg_scrcoulomb_2d_hermitian(this)
        class(idiel_t), target, intent(inout) :: this
    end subroutine compute_anisotropic_avg_scrcoulomb_2d_hermitian

end interface

end module idiel

