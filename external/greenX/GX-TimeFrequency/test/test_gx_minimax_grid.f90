! **************************************************************************************************
!  Copyright (C) 2020-2023 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************!> Note, this is essentially a unit test but there's no fortran unit test framework
!> in place. As such, we have to a) run by python, b) write to file c) python parse data and assert.
!>
!> Cannot write python bindings as all grids/weights are allocated in the library
!> - This should be changed
program test_gx_minimax_grid
    use kinds, only: dp
    use gx_minimax, only: gx_minimax_grid
    implicit none

    integer :: n_mesh_points
    real(dp) :: e_transition_min, e_transition_max
    real(dp), allocatable :: tau_mesh(:), tau_weights(:)
    real(dp), allocatable :: freq_mesh(:), freq_weights(:)
    real(dp), allocatable :: cos_tau_to_freq_weights(:, :), cos_freq_to_tau_weights(:, :)
    real(dp), allocatable :: sinft_tau_to_freq_weights(:, :)
    real(dp) :: max_errors(3)
    real(dp) :: cosft_duality_error
    integer :: ierr
    character(300) :: file_name, dummy_char, tmp_path

    integer :: i, unit_id, n

    ! Pass file name from command line
    call get_command_argument(1, file_name)

    ! Read settings from file
    open(newunit=unit_id, file=trim(adjustl(file_name)))
    ! Rather than parse through file_name to get split into tmp_path/file
    ! just put the path in the input file and extract
    read(unit_id, *) tmp_path
    read(unit_id, *) dummy_char, n_mesh_points
    read(unit_id, *) dummy_char, e_transition_min
    read(unit_id, *) dummy_char, e_transition_max
    close(unit_id)

    ! Call the library API
    call gx_minimax_grid(n_mesh_points, e_transition_min, e_transition_max, &
                         tau_mesh, tau_weights, &
                         freq_mesh, freq_weights, &
                         cos_tau_to_freq_weights, cos_freq_to_tau_weights, &
                         sinft_tau_to_freq_weights, &
                         max_errors, cosft_duality_error, ierr)

    ! Extract tmp_path from file_name, where file_name = tmp_path / 'inputs.dat'
    n = index(trim(file_name), 'inputs.dat')
    tmp_path = file_name(1:n-1)

    ! Tau grid and weights
    open(newunit=unit_id, file = trim(adjustl(tmp_path))//'tau.dat')
    do i = 1, size(tau_mesh)
        write(unit_id, *) tau_mesh(i), tau_weights(i)
    end do
    close(unit_id)

    ! Frequency grid and weights
    open(newunit=unit_id, file = trim(adjustl(tmp_path))//'freq.dat')
    do i = 1, size(freq_mesh)
        write(unit_id, *) freq_mesh(i), freq_weights(i)
    end do
    close(unit_id)

    ! TODO Finish
    ! cos and sin transform weights
!    open(newunit=unit_id, file = trim(adjustl(tmp_path))//'cos_sin_weights.dat')
!    do i = 1, size(cos_tau_to_freq_weights)
!        write(unit_id, *) cos_tau_to_freq_weights(i, j), cos_freq_to_tau_weights(i, j), &
!                    sinft_tau_to_freq_weights(i, j)
!    end do
!    close(unit_id)

    ! TODO Add max_errors and cosft_duality_error as outputs

    ! Note, arrays are allocated in the library, and so must be explicitly
    ! deallocated by the caller
    deallocate(tau_mesh, tau_weights, freq_mesh, freq_weights, &
               cos_tau_to_freq_weights, cos_freq_to_tau_weights, sinft_tau_to_freq_weights)

end program test_gx_minimax_grid
