!> This module include procedures to interpolate function $f1_n( k )$ defined on the kmesh 1, to kmesh 2
!> using 3D Smooth Fourier transform according to
!> PRB 38, 2721 (1988).
module mod_fourier_interpolation

    use precision, only: i32, dp
    use constants, only: zzero, twopi

    implicit none

    private
    public  :: smooth_fourier_interpolator

    !> Derived type holding the Fourier interpolator
    type :: smooth_fourier_interpolator
        !> Number of bands in the reference data
        integer(i32), private :: band_range(2)
        !> Number of interpolation nodes
        integer(i32), private :: n_nodes
        !> Interpolation nodes
        real(dp), allocatable, private :: interpolation_nodes(:,:)
        !> Reference data
        real(dp), allocatable, private :: reference_data(:,:)
        !> Interpolation coefficients
        complex(dp), allocatable, private :: interpolation_coefficients(:,:)
        !> Parameter C1 of the suppression function; see 10.1016/j.cpc.2006.03.007
        real(dp), private :: c1 = 0.75_dp
        !> Parameter C2 of the suppression function; see 10.1016/j.cpc.2006.03.007
        real(dp), private :: c2 = 0.75_dp
        !> Rolerance value (below which we consider zero)
        real(dp), private :: zero_tolerance = 1e-10_dp
        !> The size of the point group
        integer(i32), private :: n_point_group
    contains
        procedure, public :: initialize  => initialize_smooth_fourier_interpolator_t
        procedure, public :: clean       => clean_smooth_fourier_interpolator_t
        procedure, public :: interpolate => interpolate_smooth_fourier_interpolator_t
        procedure, public :: quality     => test_smooth_fourier_interpolator_t
    end type smooth_fourier_interpolator

contains

    !> This procedure initializes the smooth_fourier_interpolator
    subroutine initialize_smooth_fourier_interpolator_t(this, interpolation_nodes, reference_data, idx_band_start, idx_band_end)

        use mod_symmetry,      only: nsymcrys, vtlsymc
        use mod_fourier_interpolation_domain, only: init_real_space_domain_for_fourier_interpolation, &
                                                    nrr, rindex, rbas, Sigma_rindex


        !> smooth_fourier_interpolator obj to initialize
        class(smooth_fourier_interpolator), intent(out) :: this
        !> the interpolation nodes (3,npoints)
        real(dp), allocatable, intent(in)               :: interpolation_nodes(:,:)
        !> the reference data (npoints,nb)
        real(dp), allocatable, intent(in)               :: reference_data(:,:)
        !> the first band for the interpoland generation
        integer(i32), intent(in)                        :: idx_band_start
        !> the last band for the interpoland generation
        integer(i32), intent(in)                        :: idx_band_end
        
        ! Suppression function and its inverse
        real(dp), allocatable :: rho(:), rho_inv(:)

        ! Factors for the supperssion function
        real(dp) :: x2, x6

        ! R vector related variables
        real(dp) :: rlen, rlen_min, rdiv

        ! Sm matrix
        complex(dp), allocatable :: Sm(:,:), Sm_diff(:,:)

        ! Dummy indexes
        integer(i32) :: ir, inode, inodep, ipoint_group

        ! H matrix
        complex(dp), allocatable :: H(:,:)
       
        ! Energy differences matrix
        complex(dp), allocatable :: delta_e(:,:)

        ! LAPACK linear system
        integer(i32) :: info
        integer(i32), allocatable :: ipiv(:)

        ! LAPACK routine for solving a linear system
        external :: zgetrf, zgetrs

        ! Initialize variables
        this%band_range(:) = [idx_band_start, idx_band_end]
        this%n_nodes = size(interpolation_nodes, dim=2)
        allocate(this%reference_data, source=reference_data(1:this%n_nodes,idx_band_start:idx_band_end))
        allocate(this%interpolation_nodes, source=interpolation_nodes)

        ! Initialize the real space vectors for the interpoland
        call init_real_space_domain_for_fourier_interpolation()

        ! Get the number of point group symmetries (i.e. hose operations that involve only rotational components without any translational part)
        this%n_point_group = count(norm2(vtlsymc(:,1:nsymcrys),dim=1) < this%zero_tolerance)

        !----------Calculate the curvature function (rho) for each R---------------
        allocate(rho(nrr), source = 0.0_dp)
        rho(1) = 1.0_dp

        ! Get R_min. It is the second element of rindex, the first one is the null R
        rlen_min = norm2(matmul(transpose(rbas), rindex(1:3,2)))

        ! Iterate over R
        do ir = 2, nrr
            rlen = norm2(matmul(transpose(rbas), rindex(1:3,ir)))
            rdiv = rlen/rlen_min
            ! The rdiv**6 helps to suppress small amplitude, short-wavelength wiggles.
            rho(ir) = (1.0_dp - this%c2 * rdiv**2)**2 + this%c2 * rdiv**6
        end do

        !----------Compute the expansion (S_m)---------------
        allocate(Sm(this%n_nodes, nrr), source = zzero)

        !$omp parallel do collapse(2) default(none) &
        !$omp shared(this, nrr, Sm, Sigma_rindex, rindex) private(inode, ir, ipoint_group)
        do inode = 1, this%n_nodes
            do ir = 2, nrr
                do ipoint_group = 1, this%n_point_group
                    Sm(inode, ir) = Sm(inode, ir) + exp(cmplx(0.0, twopi * &
                                        dot_product(this%interpolation_nodes(1:3, inode), Sigma_rindex(1:3, ipoint_group, ir)), kind=dp))
                end do
                Sm(inode, ir) = Sm(inode, ir) / this%n_point_group
            end do
        end do
        !$omp end parallel do

        ! Clean the expression from small values and divide it by the suppression function
        Sm = merge(cmplx(0.0_dp,aimag(Sm),kind=dp), Sm, abs(real(Sm))  < this%zero_tolerance)
        Sm = merge(cmplx(real(Sm),0.0_dp,kind=dp) , Sm, abs(aimag(Sm)) < this%zero_tolerance)

        ! Now substract the last point from it
        ! Moreover compute the inverse of rho
        call move_alloc(rho, rho_inv)
        allocate(Sm_diff, source = Sm)
        do ir = 2, nrr
                Sm_diff(1:this%n_nodes-1, ir) = Sm_diff(1:this%n_nodes-1, ir) - Sm_diff(this%n_nodes, ir)
                Sm_diff(this%n_nodes, ir)     = zzero
                rho_inv(ir) = 1.0_dp / rho_inv(ir)
        end do

        ! Compute the H matrix and de delta_e vector
        allocate(H(this%n_nodes-1, nrr), source = zzero)
        allocate(delta_e(this%n_nodes-1, idx_band_start:idx_band_end), source = zzero)

        delta_e(:,:) = this%reference_data(1:this%n_nodes-1,:)

        !$omp parallel do collapse(2) default(none) &
        !$omp shared(this, nrr, H, Sm_diff, rho_inv) private(inode, inodep, ir)
        do inode = 1, this%n_nodes - 1
            do inodep = 1, this%n_nodes - 1
                do ir = 2, nrr
                    H(inode, inodep) = H(inode, inodep) + Sm_diff(inode, ir) * conjg(Sm_diff(inodep, ir)) * rho_inv(ir)
                end do
            end do
        end do
        !$omp end parallel do

        ! Solve the Linear equations for the Lagrange multipliers
        allocate(ipiv(this%n_nodes-1))
        call zgetrf(this%n_nodes - 1, this%n_nodes - 1, H, this%n_nodes - 1, ipiv, info)
        call errmsg(info .ne. 0,"initialize_smooth_fourier_interpolator_t","Error when calling zgetrf")

        call zgetrs('n', this%n_nodes - 1, this%band_range(2) - this%band_range(1) + 1, &
                     H, this%n_nodes - 1, ipiv, delta_e, this%n_nodes - 1, info)
        call errmsg(info .ne. 0,"initialize_smooth_fourier_interpolator_t","Error when calling zgetrs")

        allocate(this%interpolation_coefficients(nrr, idx_band_start:idx_band_end), source = zzero)

        ! Construct the coefficients matrix
        this%interpolation_coefficients(1,:) = this%reference_data(this%n_nodes,:)
        do ir = 2, nrr
            do inode = 1, this%n_nodes - 1
                this%interpolation_coefficients(ir, :) = this%interpolation_coefficients(ir, :) + &
                                                         delta_e(inode, :) * conjg(Sm_diff(inode, ir)) * rho_inv(ir)
            end do
            this%interpolation_coefficients(1, :) = this%interpolation_coefficients(1, :) - this%interpolation_coefficients(ir, :) * Sm(this%n_nodes, ir)
        end do

    end subroutine initialize_smooth_fourier_interpolator_t

    !> Cleans the smooth_fourier_interpolator
    subroutine clean_smooth_fourier_interpolator_t(this)
        !> smooth_fourier_interpolator obj to destroy
        class(smooth_fourier_interpolator), intent(inout) :: this
        deallocate(this%reference_data, this%interpolation_nodes, this%interpolation_coefficients)
    end subroutine clean_smooth_fourier_interpolator_t

    subroutine interpolate_smooth_fourier_interpolator_t(this, interpolation_points, interpolated_data, idx_band_start, idx_band_end)

        use mod_fourier_interpolation_domain, only: nrr, rbas, Sigma_rindex
        use general_matrix_multiplication, only: matrix_multiply

        !> smooth_fourier_interpolator obj to use
        class(smooth_fourier_interpolator), intent(in) :: this
        !> Points in which to compute the interpolation (3,npoints)
        real(dp), allocatable, intent(in) :: interpolation_points(:,:)
        !> Interpolated data
        real(dp), allocatable, intent(out) :: interpolated_data(:,:)
        !> the first band for the interpolation
        integer(i32), intent(in)                        :: idx_band_start
        !> the last band for the interpolation
        integer(i32), intent(in)                        :: idx_band_end

        integer(i32) :: npoints_interpolation, ipoint, ir, ipoint_group
        complex(dp), allocatable  :: Sm(:,:), interpolator_result(:,:)

        call errmsg(idx_band_start < this%band_range(1) .or. idx_band_start > this%band_range(2), "interpolate_smooth_fourier_interpolator_t", &
                    "The first band for the interpolation is not in the interpolator range")
        call errmsg(idx_band_end < this%band_range(1) .or. idx_band_end > this%band_range(2), "interpolate_smooth_fourier_interpolator_t", &
                    "The last band for the interpolation is not in the interpolator range")


        ! Get the number of interpolation points
        npoints_interpolation = size(interpolation_points,2)

        allocate(Sm(npoints_interpolation, nrr), source=zzero)

        !$omp parallel do collapse(2) default(none) &
        !$omp shared(npoints_interpolation, interpolation_points, this, nrr, Sm, Sigma_rindex) private(ipoint, ir, ipoint_group)
        do ipoint = 1, npoints_interpolation
            do ir = 2, nrr
                do ipoint_group = 1, this%n_point_group
                    Sm(ipoint, ir) = Sm(ipoint, ir) + exp(cmplx(0.0, twopi * &
                                        dot_product(interpolation_points(1:3, ipoint), Sigma_rindex(1:3, ipoint_group, ir)), kind=dp))
                end do
                Sm(ipoint, ir) = Sm(ipoint, ir) / this%n_point_group
            end do
        end do
        !$omp end parallel do

        allocate(interpolator_result(npoints_interpolation, idx_band_start:idx_band_end))
        call matrix_multiply(Sm, this%interpolation_coefficients, interpolator_result)

        ! Get the real part
        allocate(interpolated_data(npoints_interpolation, idx_band_start:idx_band_end))
        interpolated_data(:,:) = real(interpolator_result(:,:))

    end subroutine interpolate_smooth_fourier_interpolator_t

    !> It provides a coherence test, giving the relative difference between the reference_data
    !> and the interpolated data in the interpolation nodes.
    function test_smooth_fourier_interpolator_t(this) result(relative_difference)
        class(smooth_fourier_interpolator), intent(in) :: this

        real(dp), allocatable  :: interpolator_result(:,:)
        real(dp) :: relative_difference

        call this%interpolate(this%interpolation_nodes, interpolator_result, this%band_range(1), this%band_range(2))

        ! Get the relative difference
        relative_difference = norm2( interpolator_result - this%reference_data )   / norm2(this%reference_data)

    end function test_smooth_fourier_interpolator_t


end module mod_fourier_interpolation
