module space_time
    use precision, only: wp, dp
    use modmpi, only: mpiinfo
    use groundstate_results, only: groundstate_results_type
    use time_freq_grid, only: grid_type
    use greens_bands, only: g0_bands_img_time
    !use gx_minimax, only: gx_minimax_grid
    implicit none
    private

    public :: space_time_init, space_time_main

contains

    !> Initialise quantities required for the real-space img time algorithm
    !> from all the global exciting data.
    subroutine space_time_init(comm, gs_results, img_time_grid)
        !> MPI communicator
        type(mpiinfo), intent(inout) :: comm
        !> Ground state results
        type(groundstate_results_type), intent(out) :: gs_results
        !> Imaginary time grid
        type(grid_type), intent(out) :: img_time_grid

        !> Place holders for minimax grids/weights
        real(dp), allocatable :: dummy_points(:), dummy_weights(:)
        !>Local variables
        integer :: num_points
        real(dp) :: e_min, e_max
        real(dp), allocatable :: tau_points, tau_weights
        real(dp), allocatable :: omega_points, omega_weights
        real(dp), allocatable :: cosft_wt
        real(dp), allocatable :: cosft_tw
        real(dp), allocatable :: sinft_wt
        real(dp) :: max_errors(3)
        real(dp) :: cosft_duality_error
        integer :: ierr

        ! exciting ground state data
        call gs_results%init_from_globals()
        
        ! Should call minimax here
        allocate(dummy_points(3), source=[-1._dp, 0._dp, 1._dp])
        allocate(dummy_weights(3), source=[1._dp, 1._dp, 1._dp])
        call img_time_grid%init(dummy_points, dummy_weights, 'imaginary', label='minimax')
        ! call  gx_minimax_grid(num_points, e_min, e_max, &
        !                       tau_points, tau_weights,  & 
        !                       omega_points, omega_weights, &
        !                       cosft_wt, cosft_tw, sinft_wt, & 
        !                       max_errors, cosft_duality_error, ierr)
    end subroutine


    !> Main for space time construction of the dielectric function
    subroutine space_time_main(comm, gs_results, img_time_grid)

        !> MPI communicator
        type(mpiinfo), intent(inout) :: comm
        !> Ground state results (no idea why the compiler is complaining when I use intent(in))
        type(groundstate_results_type), intent(inout) :: gs_results
        !> Imaginary time grid
        type(grid_type), intent(in) :: img_time_grid

        !> Greens function in the band representation
        complex(wp), allocatable :: G0_bands(:, :, :)
        !> Greens function in the MT basis, real-space
        complex(wp), allocatable :: G_mt_mt_R(:, :, :, :, :)
        !> Scalar sizes
        integer :: n_ks_states, n_kpts, n_spin, n_occupied

        ! Assign scalar values from GS results object
        n_ks_states = gs_results%n_ks_states()
        n_kpts = gs_results%n_kpoints()
        n_occupied = gs_results%n_occupied()
        n_spin = gs_results%n_spin

        ! MT - MT construction of G and the polarisability 
        allocate(G0_bands(n_ks_states, n_kpts, img_time_grid%n_points()))
        call g0_bands_img_time(gs_results%eigenvalues, gs_results%efermi, img_time_grid, n_occupied, G0_bands)
        !call construct_G_realspace_coefficients(kgrid, gs_results%translations, img_time_grid, gs_results%Z, G_bands, G_mt_mt_R)
        deallocate(G0_bands)

        ! INT - INT construction of G and the polarisability 

        ! MT - INT construction of G and the polarisability 

    end subroutine

end module
