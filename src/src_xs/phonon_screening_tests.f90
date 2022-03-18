module phonon_screening_tests

    use unit_test_framework, only: unit_test_type
    use precision, only: sp, dp
    use math_utils, only: all_zero, all_close, diag
    use constants, only: zone, zzero
    use phonon_screening, only: gen_phonon_screening, select_bse_energies
    use modmpi, only: mpiinfo

    implicit none

    private

    public :: phonon_screening_test_driver

contains

    !> Runs all tests in this module
    !> Called from top-level routine
    subroutine phonon_screening_test_driver(mpiglobal, kill_on_failure)

        !> mpi information
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program before the test driver finishes,
        !> if an assertion fails
        logical, optional :: kill_on_failure

        !> Our test object that looks like the Zofu object
        type(unit_test_type) :: test_report
        !> Number of tests
        integer, parameter :: n_assertions = 9

        ! Initialise tests
        call test_report%init(n_assertions, mpiglobal)

        ! Call tests
        call test_real_diagonal(test_report)

        call test_isotropic_cubic_system(test_report)

        call test_nonpolar_system(test_report)

        call test_select_bse_energies(test_report)

        ! report results
        if (present(kill_on_failure)) then
            call test_report%report('phonon_screening', kill_on_failure)
        else
            call test_report%report('phonon_screening')
        end if

        ! Deallocate after tests
        call test_report%finalise()
    end subroutine phonon_screening_test_driver

    !> Test if the diagonal elemements of the phonon
    !> screened Coulomb interaction
    !> \[
    !>    W_{\mathbf{q}\nu}^{\mathbf{G}\mathbf{G}'}(\mathbf{q})
    !> \]
    !> are real
    subroutine test_real_diagonal(test_report)

        !> Our test object
        type(unit_test_type), intent(inout) :: test_report

        !> Tolerance of deviation for tests
        real(dp), parameter :: tol = 1.e-10_dp
        !> Fourier coeff scr. Coul. int. W_{GG'} for given phonon mode and qpoint
        complex(dp), allocatable    ::  W_GGp(:, :)
        !> Define which phonon mode should be tested
        integer(sp) :: test_mode
        !> Define which q-point should be tested
        integer(sp), parameter :: test_qpt = 1
        !> Number of phonon modes
        integer(sp) :: n_modes
        !> Total number of atoms
        integer(sp) :: n_atom_tot

        !> number of species
        integer(sp) :: n_species
        !> Number of atoms per species
        integer(sp), allocatable :: n_atom_species(:)
        !> Atomic positions
        real(dp), allocatable :: atom_positions(:, :, :)
        !> Species masses
        real(dp), allocatable :: mass_species(:)
        !> Unit cell volume
        real(dp) :: omega

        !> Number of reduced q-points
        integer(sp) :: n_qpt_r
        !> Number of k-points
        integer(sp) :: n_kpt
        !> (q+G)-vectors
        real(dp), allocatable :: q_plus_G_vec(:, :, :)
        !> Number of (q+G)-vectors
        integer(sp) :: n_q_plus_G

        !> phonon frequencies
        real(dp), allocatable    :: freq_ph(:, :)
        !> phonon eigenvectors
        complex(dp), allocatable  :: evec_ph(:, :, :, :)
        !> High frequency dielectric tensor
        real(dp) :: eps_infty(3, 3)
        !> Born-effective charge tensor
        real(dp), allocatable:: zstar(:, :, :)

        ! Set system with one atom -> 3 phonon modes
        n_species = 1
        allocate (n_atom_species(n_species), source=1)
        allocate (mass_species(n_species))
        allocate (atom_positions(3, maxval(n_atom_species),&
                                     & n_species))
        atom_positions = 1._dp
        mass_species = 1._dp
        omega = 1._dp
        n_atom_tot = sum(n_atom_species)

        ! Set 3 (q+G)-vectors:
        n_kpt = 1
        n_qpt_r = 1
        n_q_plus_G = 3
        allocate (q_plus_G_vec(3, n_q_plus_G, n_qpt_r))
        q_plus_G_vec(:, 1, :) = 1._dp
        q_plus_G_vec(:, 2, :) = 2._dp
        q_plus_G_vec(:, 3, :) = 3._dp

        ! Set  lattice screening quantities
        n_modes = 3*n_atom_tot
        allocate (freq_ph(n_modes, n_qpt_r), &
                                  & source=1._dp)
        allocate (evec_ph(3, n_atom_tot, n_modes,&
                                  &  n_qpt_r), source=cmplx(1.0, 1.0, dp))
        allocate (zstar(3, 3, n_atom_tot), source=1._dp)
        eps_infty = 1._dp

        ! Target object phonon screened Coulomb interaction
        allocate (W_GGp(n_q_plus_G, n_q_plus_G))

        ! Which mode shall be tested?
        test_mode = 1

        call gen_phonon_screening(q_plus_G_vec(:, 2, test_qpt),&
                    & q_plus_G_vec(:, :, test_qpt), &
                    & eps_infty, zstar, &
                    & freq_ph(test_mode, test_qpt),&
                    & evec_ph(:, :, test_mode, test_qpt),&
                    & n_kpt, omega, n_atom_species,&
                    & mass_species, atom_positions, W_GGp)

        ! Check if imaginary parts of all diagonal elements are 0
        call test_report%assert(all_zero(aimag(diag(W_GGp)), tol), &
                              & "Test gen_polar_phonon_screening. &
                              &  Expected: Diagonal elements are always real.")

    end subroutine test_real_diagonal

    !> Test if the phonon screened Coulomb interaction
    !> \[
    !>    W_{\mathbf{q}\nu}^{\mathbf{G}\mathbf{G}'}(\mathbf{q})
    !> \]
    !> vanishes in nonpolar materials with zero Born-effective charges
    subroutine test_nonpolar_system(test_report)

        !> Our test object
        type(unit_test_type), intent(inout) :: test_report

        !> number of species
        integer(sp) :: n_species
        !> Number of atoms per species
        integer(sp), allocatable :: n_atom_species(:)
        !> Atomic positions
        real(dp), allocatable :: atom_positions(:, :, :)
        !> Species masses
        real(dp), allocatable :: mass_species(:)
        !> Unit cell volume
        real(dp) :: omega

        !> Number of reduced q-points
        integer(sp) :: n_qpt_r
        !> Number of k-points
        integer(sp) :: n_kpt
        !> (q+G)-vectors
        real(dp), allocatable :: q_plus_G_vec(:, :, :)
        !> Number of (q+G)-vectors
        integer(sp) :: n_q_plus_G

        !> phonon frequencies
        real(dp), allocatable    :: freq_ph(:, :)
        !> phonon eigenvectors
        complex(dp), allocatable  :: evec_ph(:, :, :, :)
        !> High frequency dielectric tensor
        real(dp) :: eps_infty(3, 3)
        !> Born-effective charge tensor
        real(dp), allocatable:: zstar(:, :, :)

        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Fourier coeff scr. Coul. int. W_{GG'} for given phonon mode and qpoint
        complex(dp), allocatable    ::  W_GGp(:, :)
        !> Define which phonon mode should be tested
        integer(sp) :: test_mode
        !> Define which q-point should be tested
        integer(sp), parameter :: test_qpt = 1
        !> Number of phonon modes
        integer(sp) :: n_modes
        !> Total number of atoms
        integer(sp) :: n_atom_tot

        ! Set system with one atom -> 3 phonon modes
        n_species = 1
        allocate (n_atom_species(n_species), source=1)
        allocate (mass_species(maxval(n_atom_species)))
        allocate (atom_positions(3, maxval(n_atom_species),&
                                     & n_species))
        atom_positions = 1._dp
        mass_species = 1._dp
        omega = 1._dp
        n_atom_tot = sum(n_atom_species)
        ! Set 3 (q+G)-vectors:
        n_kpt = 1
        n_qpt_r = 1
        n_q_plus_G = 3
        allocate (q_plus_G_vec(3, n_q_plus_G, n_qpt_r))
        q_plus_G_vec(:, 1, :) = 1._dp
        q_plus_G_vec(:, 2, :) = 2._dp
        q_plus_G_vec(:, 3, :) = 3._dp

        ! Set all  phonon quantities to 1
        n_modes = 3*n_atom_tot
        allocate (freq_ph(n_modes, n_qpt_r),&
                                  & source=1._dp)
        allocate (evec_ph(3, n_atom_tot, n_modes,&
                                    &n_qpt_r), source=cmplx(1.0, 1.0, dp))
        allocate (zstar(3, 3, n_atom_tot), source=0._dp)
        eps_infty = 1._dp

        ! Target object phonon screened Coulomb interaction
        allocate (W_GGp(n_q_plus_G, n_q_plus_G))

        ! Which mode shall be tested?
        test_mode = 1
        call gen_phonon_screening(q_plus_G_vec(:, 2, test_qpt),&
                &q_plus_G_vec(:, :, test_qpt),&
                & eps_infty, zstar, &
                & freq_ph(test_mode, test_qpt),&
                & evec_ph(:, :, test_mode, test_qpt),&
                & n_kpt, omega, n_atom_species,&
                & mass_species, atom_positions, W_GGp)

        ! Check screened Coulomb vanishes
        call test_report%assert(all_zero(W_GGp, tol), &
                              & "Test gen_polar_phonon_screening. &
                              &  Expected: Zero contribution in nonpolar system")

    end subroutine test_nonpolar_system

    !> 3 Tests for a cubic and isotropic system in a sense that
    !> the Born-effecticve charge tensor \[ \mathbf{Z}^*_\kappa \]
    !> and the dielectric tensor \[ \boldsymbol{\varepsilon}_\infty \]
    !> are diagonal having only one independent element.
    !> Checks that the following conditions are fullfilled:
    !> 1) Zero contribution if phonon eigenvector is perpendicular to (q+G)
    !> 2) Non-zero contribution if phonon eigenvector is parallel to (q+G)
    !> 3) W_{GG'}(q)*|G+q||G'+q| is independent of G and G' for one q
    !>    if phonon eigenvector is parallel to (q+G) and (q+G')
    !>    i.e. W_{GG'}(q) should be of the form W_{GG'}(q) = f(q)*1/|G+q|*1/|G'+q|
    subroutine test_isotropic_cubic_system(test_report)

        !> Our test object
        type(unit_test_type), intent(inout) :: test_report

        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Fourier coeff scr. Coul. int. W_{GG'} for given phonon mode and qpoint
        complex(dp), allocatable    ::  W_GGp(:, :)
        !> Define index of longitudinal mode
        integer(sp), parameter :: long_mode = 3
        !> Define index transveres mode
        integer(sp), parameter :: trans_mode = 2
        !> Running indeces
        integer(sp) :: i, j
        !> Matrix filled with ones
        complex(dp), allocatable :: one_matrix(:, :)
        !> Help matrix for tests
        complex(dp), allocatable :: aux_matrix(:, :)
        !> Define which phonon mode should be tested
        integer(sp) :: test_mode
        !> Define which q-point should be tested
        integer(sp), parameter :: test_qpt = 1
        !> Number of phonon modes
        integer(sp) :: n_modes
        !> Total number of atoms
        integer(sp) :: n_atom_tot

        !> number of species
        integer(sp) :: n_species
        !> Number of atoms per species
        integer(sp), allocatable :: n_atom_species(:)
        !> Atomic positions
        real(dp), allocatable :: atom_positions(:, :, :)
        !> Species masses
        real(dp), allocatable :: mass_species(:)
        !> Unit cell volume
        real(dp) :: omega

        !> Number of reduced q-points
        integer(sp) :: n_qpt_r
        !> Number of k-points
        integer(sp) :: n_kpt
        !> (q+G)-vectors
        real(dp), allocatable :: q_plus_G_vec(:, :, :)
        !> Number of (q+G)-vectors
        integer(sp) :: n_q_plus_G

        !> phonon frequencies
        real(dp), allocatable    :: freq_ph(:, :)
        !> phonon eigenvectors
        complex(dp), allocatable  :: evec_ph(:, :, :, :)
        !> High frequency dielectric tensor
        real(dp) :: eps_infty(3, 3)
        !> Born-effective charge tensor
        real(dp), allocatable:: zstar(:, :, :)

        ! Define model system with one atom -> 3 phonon modes
        n_species = 1
        allocate (n_atom_species(1), source=1)
        allocate (mass_species(maxval(n_atom_species)))
        allocate (atom_positions(3, maxval(n_atom_species),&
                                     & n_species))

        ! Set system with atoms in the origin
        atom_positions = 0._dp
        mass_species = 1._dp
        omega = 1._dp

        ! Set 3 (q+G)-vectors with components only in x-direction:
        n_kpt = 1
        n_qpt_r = 1
        n_q_plus_G = 3
        allocate (q_plus_G_vec(3, n_q_plus_G, n_qpt_r), source=0._dp)
        q_plus_G_vec(1, :, test_qpt) = [1._dp, 2._dp, 3._dp]

        ! Set longitudinal mode non-zero in x-direction (i.e. parallel to all (q+G))
        ! Set transverse mode non-zero in y-direction (perpendicular to all (q+G))
        n_atom_tot = sum(n_atom_species)
        n_modes = 3*n_atom_tot
        allocate (freq_ph(n_modes, n_qpt_r))
        allocate (evec_ph(3, n_atom_tot, n_modes,&
                                   & n_qpt_r), source=zzero)
        freq_ph = 1._dp
        evec_ph(1, :, long_mode, :) = 1._dp
        evec_ph(2, :, trans_mode, :) = 1._dp

        ! Use diagonal eps_infty and Born-effective charges charges
        allocate (zstar(3, 3, n_atom_tot), source=0._dp)
        eps_infty = 0._dp
        do j = 1, 3
            eps_infty(j, j) = j*1._dp
            zstar(j, j, 1) = j*1._dp
        end do

        ! Target object W_GGp and auxiliary matrices
        allocate (W_GGp(n_q_plus_G, n_q_plus_G))
        allocate (aux_matrix(n_q_plus_G, n_q_plus_G))
        allocate (one_matrix(n_q_plus_G, n_q_plus_G), source=zone)
        evec_ph(1, 1, long_mode, 1) = 1._dp
        evec_ph(2, 1, trans_mode, 1) = 1._dp

        ! TEST 1
        call gen_phonon_screening(q_plus_G_vec(:, 2, test_qpt),&
                & q_plus_G_vec(:, :, test_qpt), eps_infty,&
                & zstar, freq_ph(trans_mode, test_qpt),&
                & evec_ph(:, :, trans_mode, test_qpt), n_kpt, &
                & omega, n_atom_species, mass_species,&
                & atom_positions, W_GGp)

        call test_report%assert(all_zero(W_GGp, tol), &
                              & "Test gen_polar_phonon_screening. &
                              &  Expected: Transverse mode does not contribute")

        ! TEST 2
        call gen_phonon_screening(q_plus_G_vec(:, 2, test_qpt),&
                & q_plus_G_vec(:, :, test_qpt), eps_infty,&
                & zstar, freq_ph(long_mode, test_qpt),&
                & evec_ph(:, :, long_mode, test_qpt),&
                & n_kpt, omega, n_atom_species,&
                & mass_species, atom_positions, W_GGp)

        call test_report%assert(.not. all_zero(W_GGp, tol), &
                              & "Test gen_polar_phonon_screening. &
                              &  Expected: Longitudinal mode contributes")

        ! TEST 3 (re-uses output of TEST 2)
        do i = 1, n_q_plus_G
            do j = 1, n_q_plus_G
                aux_matrix(i, j) = W_GGp(i, j)*norm2(q_plus_G_vec(:, i, test_qpt))&
                                  &*norm2(q_plus_G_vec(:, j, test_qpt))
            end do
        end do

        ! Check if all elements of auxiliary matrix are equal
        call test_report%assert(all_close(aux_matrix, &
                            & aux_matrix(1, 1)*one_matrix, tol), &
                            & "Test gen_polar_phonon_screening. &
                            &  Expected: Screened Coulomb of form &
                            &  W_{GG'}(q) = f(q)*1/|G+q|*1/|G'+q|")

    end subroutine test_isotropic_cubic_system

    !> Tests the selection of KS eigenenergies for BSE
    subroutine test_select_bse_energies(test_report)

        !> Our test object
        type(unit_test_type), intent(inout) :: test_report
        !> Upper/lower bounds for valence and conduction bands
        integer(sp) :: limits_bse(4)
        !> Number of unoccupied bands considered in BSE
        integer(sp) :: n_unocc_bands
        !> Number of occupied bands considered in BSE
        integer(sp) :: n_occ_bands
        !> Total number computed KS bands
        integer(sp) :: n_ks_bands
        !> Scissor paratemer
        real(dp) :: scissor
        !> 10 computed Kohn Sham energy eigenvalues at k'
        real(dp) :: eigvals_k(6)
        !> 6 computed Kohn Sham energy eigenvalues at k'
        real(dp) :: eigvals_kp(6)
        !> Scissor parameter
        real(dp) :: bse_scissor
        !> Energy eigenvalues of occupied the states considered in BSE at k,k'
        real(dp), allocatable :: eigvals_occ_bse(:, :)
        !> Energy eigenvalues of unoccupied states considered in BSE at k,k'
        real(dp), allocatable :: eigvals_unocc_bse(:, :)
        !> Running index
        integer(sp) :: i
        !> Store results
        logical :: results(4)
        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Lowest occupied state (absolut index)
        integer(sp) ::  lowest_occ_state
        !> Highest occupied state (absolut index)
        integer(sp) ::  highest_occ_state
        !> Lowest occupied state (absolut index)
        integer(sp) ::  lowest_unocc_state
        !> Lowest occupied state (absolut index)
        integer(sp) ::  highest_unocc_state

        lowest_occ_state = 2
        highest_occ_state = 3
        lowest_unocc_state = 4
        highest_unocc_state = 5
        limits_bse = [lowest_unocc_state, highest_unocc_state,&
                    & lowest_occ_state, highest_occ_state]
        n_unocc_bands = highest_unocc_state - lowest_unocc_state + 1
        n_occ_bands = highest_occ_state - lowest_occ_state + 1
        n_ks_bands = size(eigvals_kp)

        allocate (eigvals_occ_bse(n_occ_bands, 2))
        allocate (eigvals_unocc_bse(n_unocc_bands, 2))

        !Fill KS eigenvalues
        do i = 1, n_ks_bands
            eigvals_k(i) = i*1._dp
            eigvals_kp(i) = i*1._dp + 0.1_dp
        end do
        bse_scissor = 0.5_dp

        call select_bse_energies(limits_bse, eigvals_k, eigvals_kp, bse_scissor,&
                          &  eigvals_occ_bse, eigvals_unocc_bse)

        call test_report%assert(all_close(eigvals_occ_bse(:, 1),&
                & [2._dp, 3._dp], tol), &
                & 'Test select_bse_energies. Expected: [2._dp,3._dp]')
        call test_report%assert(all_close(eigvals_occ_bse(:, 2),&
                & [2.1_dp, 3.1_dp], tol),&
                & 'Test select_bse_energies. Expected: [2.1_dp,3.1_dp]')
        call test_report%assert(all_close(eigvals_unocc_bse(:, 1), &
                & [4.5_dp, 5.5_dp], tol),&
                & 'Test select_bse_energies. Expected: [4.5_dp,5.5_dp]')
        call test_report%assert(all_close(eigvals_unocc_bse(:, 2),&
                & [4.6_dp, 5.6_dp], tol),&
                & 'Test select_bse_energies. Expected: [4.6_dp,5.6_dp]')

    end subroutine test_select_bse_energies

end module phonon_screening_tests
