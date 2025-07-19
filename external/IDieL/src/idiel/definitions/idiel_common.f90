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
submodule (idiel) idiel_common

    use iso_c_binding
#include "offload.fpp"

    implicit none

contains

    module subroutine init_common(this, lattice, redpos, elements, nq, dim, nsym, crot, device_world, host_world)

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

        ! Locals
        integer(i32) :: ii, ll, mm
        real(aip) :: v_bz
        real(aip) :: rel_error
        real(aip), allocatable :: xyz(:,:)

        ! The volume in which the integral is performed
        real(aip) :: v_integral
        ! The distance to subcell surface with Gamma in the center
        real(aip), allocatable :: kmax(:)
        ! Geometric part of the integral at given point, note that is multiplied by the weight
        real(aip), allocatable :: kmax_f(:)
        ! Spherical harmonics
        complex(aip), allocatable :: ylm(:,:), ylm_pair(:,:)
        ! Some cutoff values to create the radial mesh for the 2D case
        real(aip) :: rmax, dx
        ! The reciprocal vectors
        real(aip) :: a(3), b(3)

        ! Initialize the crystal structure
        call this%cell%initialize(lattice, redpos, elements)

        ! Init reciprocal mesh size
        this%nq = nq

        ! Initialize the symmetry using the external code if not fallback to SPGLIB (if available otherwise we raise error)
        if (present(nsym) .and. present(crot)) then
            this%symmetry%nsym = nsym
            allocate(this%symmetry%crot(3,3,nsym))
            this%symmetry%crot = crot
        else
#ifdef USE_SPGLIB
            call this%symmetry%initialize(this%cell)      
#else 
            error stop "Error(idiel_t%init_common): missing symmetry information"
#endif  
        endif

        ! Process the GPU-CPU communicator
        if (present(device_world)) then
            this%world => device_world
            this%internal_world = .false.
        else
#ifdef DEVICEOFFLOAD
            if (.not. present(host_world)) error stop "Error(idiel_t%init_common): For device either the MPI communicator or the device communicator are needed"
            allocate(this%world)
            call this%world%init(host_world)
#else
            allocate(this%world)
            call this%world%init(0)
#endif
            this%internal_world = .true.
        end if


        ! Process dimensionality
        if (present(dim)) this%dim = dim

        select case(this%dim)
        
        case(2) ! 2D case

            ! rcut
            this%rcut = 0.5_aip * this%cell%lattice(3,3) 

            ! Initalize the big angular mesh and the weights for the angular part of the integral
            call compute_angular_mesh_gauss_legendre(size_mesh_2d_fine , this%ang, this%weights_fine, xyz)

            ! Add the area factor to the weights 
            a = this%cell%rlattice(1,:) / this%nq(1)
            b = this%cell%rlattice(2,:) / this%nq(2)
            this%weights_fine = this%weights_fine / & 
                sqrt( (a(2)*b(3) - a(3)*b(2))**2 +  (a(3)*b(1) - a(1)*b(3))**2 + (a(1)*b(2) - a(2)*b(1))**2 )

            ! Save the large phi mesh
            allocate(this%phi, source=this%ang(:,2))
            
            ! Compute circular basis in the fine mesh
            call circ_harm(lmax, this%ang(:,2), this%blm_fine)

            ! Get the number of points used in the angular quadrature
            this%quadrature_npoints = size(xyz, 1)

            ! Compute kmax (the boundary)
            call this%cell%get_kmax_subcell_bz(this%nq, xyz, this%rmax2d)

            ! Init the radial mesh (this mesh is intended to provide good treatment of the queue)
            rmax = 1.02_aip * maxval(this%rmax2d)
            dx = rmax / (nr - 1)
            allocate(this%radii(nr))
            do ii = 1, nr
                this%radii(ii) = (ii - 1) * dx
            end do

            ! Precompte exponential factors
            allocate(this%vr(size_mesh_2d_coarse,nr))
            do ii = 1, nr
                this%vr(:,ii) = fourpi * (1.0_aip - exp(-this%rcut * this%radii(ii)))
            end do

            ! Compute the small mesh
            deallocate(this%ang, xyz)
            call compute_angular_mesh_gauss_legendre(size_mesh_2d_coarse, this%ang, this%weights, xyz)
            allocate(this%xyz, source=transpose(cmplx(xyz,0.0,aip)))
            this%quadrature_npoints = size(xyz, 1)
            ! Compute circular basis in the coarse mesh
            call circ_harm(lmax, this%ang(:,2), this%blm_coarse)

            OMP_OFFLOAD target enter data map(to: this%xyz, this%weights, this%blm_fine, this%weights_fine, this%radii, this%rmax2d, this%vr)

        case(3) ! 3D case
            ! Compute the volume in which the integral will be performed
            v_bz = 8.0 * pi**3 / this%cell%vuc 
            v_integral = product(this%nq) / v_bz

            ! Initalize the big angular mesh and the weights for the angular part of the integral
            call compute_angular_mesh_lebedev_131(this%ang, this%weights_fine, xyz)

            ! Get the number of points used in the angular quadrature
            this%quadrature_npoints = size(xyz, 1)

            ! Compute kmax (the boundary)
            call this%cell%get_kmax_subcell_bz(this%nq, xyz, kmax)

            ! Construct the geometric function
            allocate(kmax_f(this%quadrature_npoints))

            ! Compute the geometric part of the integral
            ! Notice that this needs to be done in a dense mesh since
            ! we are trying to integrate parallelepiped using spherical coordinates 
            ! and that is an object of infinite degree as for its spherical harmonic expansion
            ! K function :  \frac{q_{max}^{3}}{3V_{\Gamma}} 
            kmax_f(:) =  this%weights_fine * v_integral * onethird * kmax(:)**3

            ! Precompute the angular part, we first need the spherical harmonics in the big mesh
            call sph_harm(lmax, this%ang, ylm)

            allocate(ylm_pair(size(ylm,1),nsph_pair))
            ylm_pair(:,1) = ylm(:,1)
            ii = 2
            do ll = 1, lmax
                if ( modulo(ll, 2_i32) /= 0) cycle
                do mm = -ll, ll
                    ylm_pair(:,ii) = ylm(:,ll**2 + mm + ll + 1)
                    ii = ii + 1
                end do
            end do

            allocate(this%angular_integrals(nsph_pair))

            !$omp parallel shared(ylm_pair, this, kmax_f) private(ii)
            !$omp do schedule(dynamic)
            do ii = 1, nsph_pair
                this%angular_integrals(ii) = sum(ylm_pair(:, ii) * kmax_f(:))
            end do 
            !$omp end do
            !$omp end parallel

            ! Set now to smaller mesh, as dielectric terms converge quite fast with l in comparison
            ! to the geometric part, and thus a smaller Lebedev order can be used
            deallocate(this%ang, this%weights_fine, xyz, ylm, ylm_pair)
            ! Recompute things in the small mesh
            call compute_angular_mesh_lebedev_41(this%ang, this%weights, xyz)
            allocate(this%xyz,source=transpose(cmplx(xyz,0.0,aip)))
            this%quadrature_npoints = size(xyz, 1)

            ! Compute the spherical harmonics in the smaller mesh (save only the pair values)
            call sph_harm(lmax, this%ang, ylm)

            allocate(this%ylm(size(ylm,1),nsph_pair))
            this%ylm(:,1) = ylm(:,1)
            ii = 2
            do ll = 1, lmax
                if ( modulo(ll, 2_i32) /= 0_i32) cycle
                do mm = -ll, ll
                    this%ylm(:, ii) = ylm(:,ll**2 + mm + ll + 1)
                    ii = ii + 1
                end do
            end do

            deallocate(ylm, xyz)

            OMP_OFFLOAD target enter data map(to: this%xyz, this%weights, this%ylm, this%angular_integrals)

        case default
            error stop "Error(idiel_t%init_common): Error dimension should be either 3 or 2"
        end select 
        
    end subroutine init_common

    !> This nullify and deallocates the objects
    !> @param[in] this - idiel_t object
    module subroutine clean(this)

        class(idiel_t), intent(inout) :: this

        ! Clean the device objects
        OMP_OFFLOAD target exit data map(delete: this%xyz, this%weights, this%ylm, this%angular_integrals) if(this%dim == 3)
        OMP_OFFLOAD target exit data map(delete: this%xyz, this%weights, this%blm_fine, this%weights_fine, this%radii, this%rmax2d, this%vr) if(this%dim == 2)

        if (this%internal_world) then
            call this%world%finish()
            deallocate(this%world)
        else
            nullify(this%world)
        end if

        if (associated(this%head))        nullify(this%head)
        if (associated(this%wingL))       nullify(this%wingL)
        if (associated(this%wingU))       nullify(this%wingU)
        if (associated(this%Binv))        nullify(this%Binv)
        if (allocated(this%Binv_data))    deallocate(this%Binv_data)
        if (allocated(this%weights_fine)) deallocate(this%weights_fine)
        if (allocated(this%weights))      deallocate(this%weights)
        if (allocated(this%xyz))          deallocate(this%xyz)
        if (allocated(this%ang))          deallocate(this%ang)
        if (allocated(this%idiel_wingL))  deallocate(this%idiel_wingL)
        if (allocated(this%idiel_wingU))  deallocate(this%idiel_wingU)
        if (allocated(this%idiel_body))   deallocate(this%idiel_body)
        
        ! 2D stuff
        if (allocated(this%phi))         deallocate(this%phi)
        if (allocated(this%rmax2d))      deallocate(this%rmax2d)
        if (allocated(this%blm_coarse))  deallocate(this%blm_coarse)
        if (allocated(this%blm_fine))    deallocate(this%blm_fine)
        if (allocated(this%vr))          deallocate(this%vr)
        if (allocated(this%radii))       deallocate(this%radii)
        
        ! 3D stuff
        if (allocated(this%angular_integrals)) deallocate(this%angular_integrals)
        if (allocated(this%ylm)) deallocate(this%ylm)
        
    end subroutine clean

    !> This subroutine inits pointers to a given value
    !> @param[in] this - idiel_t object
    !> @param[in] h    - head of the dielectric matrix
    !> @param[in] wl   - lower wing of the dielectric matrix
    !> @param[in] wu   - upper wing of the dielectric matrix
    !> @param[in] ib   - inverse of the body of the dielectric matrix (optional)
    module subroutine set_dielectric_blocks(this, h, wl, wu, ib)
        
        class(idiel_t), intent(inout) :: this
        complex(aip), target, intent(in)      :: h(:,:)
        complex(aip), target, intent(in)      :: wl(:,:)
        complex(aip), target, intent(in)      :: wu(:,:)
        complex(aip), target, optional, intent(in) :: ib(:,:)

        ! Associating to internal objects
        this%head    =>  h
        this%wingL   =>  wl
        this%wingU   =>  wu
        if (present(ib)) this%Binv    =>  ib

    end subroutine set_dielectric_blocks


    module subroutine invert_body(this, body)
        
        use idiel_linalg, only: inverse_complex_LU

        class(idiel_t), intent(inout), target :: this
        complex(aip), intent(in) :: body(:,:)
        
        ! Clean 
        if (associated(this%Binv)) nullify(this%Binv)
        if (allocated(this%Binv_data)) deallocate(this%Binv_data)
        
        ! Compute
        call inverse_complex_LU(body, this%Binv_data, this%world)
        
        ! Associate vector    
        this%Binv => this%Binv_data

    end subroutine invert_body

    module function get_n_basis(this) result(nbasis)
        class(idiel_t), intent(inout), target :: this
        integer(i32) :: nbasis
        if (.not. associated(this%Binv)) error stop "idiel_t%get_n_basis: Error set inverse dielectric matrix for this" 
        nbasis = size(this%Binv, 1)
    end function

end submodule idiel_common
