module putgeteps0_tests

    use unit_test_framework, only: unit_test_type
    use precision, only: sp, dp
    use math_utils, only: all_close
    use constants, only: zone, zzero
    use putgeteps0, only: geteps0_finite_q, geteps0_zero_q, &
                            puteps0_finite_q, puteps0_zero_q
    use modmpi, only: mpiinfo
    use mock_arrays, only: fill_array, value_map_real_rank1,&
                            value_map_complex_rank2,&
                            value_map_complex_rank3,&
                            value_map_complex_rank4

    implicit none

    private

    public :: putgeteps0_test_driver

contains

    !> Runs all tests in this module
    !> Called from top-level routine
    subroutine putgeteps0_test_driver(mpiglobal, kill_on_failure)

        !> mpi information
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program before the test driver finishes,
        !> if an assertion fails
        logical, optional :: kill_on_failure
        !> System command
        character(256) :: syscommand
        !> Temporary directory name
        character(256), parameter :: tmpdirname = 'tmpdir'

        !> Our test object that looks like the Zofu object
        type(unit_test_type) :: test_report
        !> Number of tests
        integer, parameter :: n_assertions = 4

        ! Initialise tests
        call test_report%init(n_assertions, mpiglobal)

        ! TODO (Alex, Max) see #127: Add tmpdir object
        ! Create temporary directory for I/O
        syscommand = 'mkdir '//trim(tmpdirname)
        call system(trim(adjustl(syscommand)))

        ! Call tests
        call test_putgeteps0_finite_q(test_report)

        call test_putgeteps0_zero_q(test_report)

        ! report results
        if (present(kill_on_failure)) then
            call test_report%report('putgeteps0', kill_on_failure)
        else
            call test_report%report('putgeteps0')
        end if

        ! Remove temporary directory
        syscommand = 'rm -r '//trim(tmpdirname)
        call system(trim(adjustl(syscommand)))

        ! Deallocate after tests
        call test_report%finalise()
    end subroutine putgeteps0_test_driver


    !> Test puteps0_finite_q and geteps0_finite_q, i.e.
    !> I/O of the dielectric matrix 
    !> \[
    !>   \tilde{\epsilon}^0_{\bf{GG'}}({\bf q},\omega) 
    !>                                                         \]
    !> for \( \mathbf{q} \neq 0\).
    subroutine test_putgeteps0_finite_q(test_report)
        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report
        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Real test input
        real(dp) :: a(3, 4)
        !> Number of G+q-vectors
        integer(sp), parameter :: n_gqvecs = 5
        !> Number of frequencies
        integer(sp), parameter :: n_freqs = 4
        !> Reference body of dielectric matrix
        complex(dp) :: eps_body(n_gqvecs, n_gqvecs, n_freqs)
        !> Body of dielectric matrix as read from file
        complex(dp) :: eps_read(n_gqvecs, n_gqvecs, n_freqs)
        !> Running index frequencies
        integer(sp) :: i_freq
        !> Frequencies
        real(dp) :: freqs(n_freqs)
        !> Filename for dielectric matrix
        character(256), parameter :: filename = 'tmpdir/EPS0_finite_q.OUT'
        !> q-vector
        real(dp):: qvec(3)

        ! Fill array with numbers according to value map
        call fill_array(eps_body, value_map_complex_rank3)
        call fill_array(qvec, value_map_real_rank1)
        call fill_array(freqs, value_map_real_rank1)

        ! Write dielectric matrix to file
        do i_freq = 1, n_freqs
            call puteps0_finite_q(1, qvec, i_freq, freqs(i_freq), &
                                  eps_body(:, :, i_freq), fname=trim(filename))
        end do

        ! Read dielectric matrix from file
        do i_freq = 1, n_freqs
            call geteps0_finite_q(1, qvec, i_freq, freqs(i_freq), &
                                  eps_read(:, :, i_freq), fname=trim(filename))
        end do

        ! Compare written and read data
        call test_report%assert(all_close(eps_body, eps_read), &
                                'Test putgeteps0_finite_q: &
                                Expected: Body of dielectric matrix as read from file equals the one that was written to file.')

    end subroutine test_putgeteps0_finite_q


    !> Test puteps0_zero_q and geteps0_zero_q, i.e.
    !> I/O of the dielectric matrix 
    !> \[
    !>   \tilde{\epsilon}^0_{\bf{GG'}}({\bf q},\omega) 
    !>                                                         \]
    !> for \( \mathbf{q} = 0\).
    !> In contrast to finite q-vectors, the dielectric matrix
    !> has head and wings.
    subroutine test_putgeteps0_zero_q(test_report)

        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report

        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Real test input
        real(dp) :: a(3, 4)
        !> Number of G+q-vectors
        integer(sp), parameter :: n_gqvecs = 5
        !> Number of frequencies
        integer(sp), parameter :: n_freqs = 4
        !> Reference body of dielectric matrix
        complex(dp) :: eps_body(n_gqvecs, n_gqvecs, n_freqs)
        !> Body of dielectric matrix as read from file
        complex(dp) :: eps_body_read(n_gqvecs, n_gqvecs, n_freqs)
        !> Reference head of dielectric matrix
        complex(dp) :: eps_head(3, 3, n_freqs)
        !> Head of dielectric matrix as read from file
        complex(dp) :: eps_head_read(3, 3, n_freqs)
        !> Reference wings of dielectric matrix
        complex(dp) :: eps_wings(n_gqvecs, 2, 3, n_freqs)
        !> Wings of dielectric matrix as read from file
        complex(dp) :: eps_wings_read(n_gqvecs, 2, 3, n_freqs)
        !> Running index frequencies
        integer(sp) :: i_freq
        !> Frequencies
        real(dp) :: freqs(n_freqs)
        !> Filename for dielectric matrix
        character(256), parameter :: filename = 'tmpdir/EPS0_zero_q.OUT'
        !> q-vector
        real(dp):: qvec(3)

        ! Fill array with numbers according to value map
        call fill_array(eps_head, value_map_complex_rank3)
        call fill_array(eps_wings, value_map_complex_rank4)
        call fill_array(eps_body, value_map_complex_rank3)
        call fill_array(qvec, value_map_real_rank1)
        call fill_array(freqs, value_map_real_rank1)

        ! Write dielectric matrix to file
        do i_freq = 1, n_freqs
            call puteps0_zero_q(1, qvec, i_freq, freqs(i_freq), &
                                eps_body(:, :, i_freq), &
                                eps_wings(:, :, :, i_freq), &
                                eps_head(:, :, i_freq), &
                               fname=trim(filename))
        end do

        ! Read dielectric matrix from file
        do i_freq = 1, n_freqs
            call geteps0_zero_q(1, qvec, i_freq, freqs(i_freq), &
                                eps_body_read(:, :, i_freq), &
                                eps_wings_read(:, :, :, i_freq), &
                                eps_head_read(:, :, i_freq), &
                                fname=trim(filename))
        end do

        ! Compare written and read data
        call test_report%assert(all_close(eps_body, eps_body_read), &
                                'Test putgeteps0_zero_q: &
                                Expected: Body of dielectric matrix as read from file &
                                equals the one that was written to file.')

        call test_report%assert(all_close(eps_head, eps_head_read), &
                                'Test putgeteps0_zero_q: &
                                Expected: Head of dielectric matrix as read from file &
                                equals the one that was written to file.')

        call test_report%assert(all_close(eps_wings, eps_wings_read), &
                                'Test putgeteps0_zero_q: &
                                Expected: Wings of dielectric matrix as read from file &
                                equals the one that was written to file.')

    end subroutine test_putgeteps0_zero_q

end module putgeteps0_tests
