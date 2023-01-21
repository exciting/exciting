module expand_add_eps_tests

    use unit_test_framework, only: unit_test_type
    use precision, only: sp, dp
    use math_utils, only: all_close
    use constants, only: zone, zzero
    use modmpi, only: mpiinfo
    use expand_add_eps_helper, only:  fold_uc_to_sc_qgrid
    use diel_mat_type, only: dielectric_matrix_type
    use expand_add_eps, only: expand_epsilon

    implicit none

    private

    public :: expand_add_eps_test_driver

contains

    !> Runs all tests in this module
    !> Called from top-level routine
    subroutine expand_add_eps_test_driver(mpiglobal, kill_on_failure)

        !> mpi information
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program before the test driver finishes,
        !> if an assertion fails
        logical, optional :: kill_on_failure

        !> Our test object that looks like the Zofu object
        type(unit_test_type) :: test_report
        !> Number of tests
        integer, parameter :: n_assertions = 3
        ! Initialise tests
        call test_report%init(n_assertions, mpiglobal)

        call test_expand_epsilon(test_report)

        call test_fold_uc_to_sc_qgrid(test_report)

        ! report results
        if (present(kill_on_failure)) then
            call test_report%report('expand_add_eps', kill_on_failure)
        else
            call test_report%report('expand_add_eps')
        end if

        ! Deallocate after tests
        call test_report%finalise()
    end subroutine expand_add_eps_test_driver

    !> Test the expansion of the dielectric matrix computed for a unit cell
    !> to the dielectric matrix of 2 x 2 supercell. Test for 2d-materials,
    !> i.e. expansion of q-vectors lying only in the x-y plane.
    !> The test is structured as follows:
    !>  i) 3 reciprocal lattice vectors of the supercell are defined.
    !>  ii) A set of 2 supercell q-vectors is defined.
    !>  iii) The set of 8 unit-cell q-vectors is constructed by shifting
    !>       the supercell q-vectors by the supercell rec. lattice vectors.
    !>  iv) The (G+q)-vectors of the supercell and of the unit-cell are
    !>      equivalent to the q-vectors of the unit-cell, the number of
    !>      (G+q)-vectors per q-vector is therefore 4  for the supercell
    !>      and 1 for the unit-cell (the q-vectors itself).
    !>  v) The diagonal unit-cell dielectric matrix is mocked according to
    !>
    !> \[
    !>       \varepsilon_{\mathbf{G}_{\rm uc}\mathbf{G}_{\rm uc}} (\mathbf{q}_{\rm uc}) =
    !>                      |\mathbf{G}_{\rm uc}|
    !>
    !>                                                                  \].
    !>  vi) The unit-cell matrix is expanded to the one of the supercell.
    !>  vii) After the expansion, it is tested if only those matrix elements
    !>       of the supercell are non-zero for which
    !>       \[
    !>          \mathbf{G}_{\rm uc} + \mathbf{q}_{\rm uc} =
    !>          \mathbf{G}_{\rm sc} + \mathbf{q}_{\rm sc}
    !>
    !>                                                           \]
    !> holds.
    subroutine test_expand_epsilon(test_report)

        use constants, only: zone
        use mod_kpointset, only: Gk_set, k_set
        use grid_utils, only: indices_zero_vectors, indices_finite_vectors

        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report
        !> Unit cell q-grid
        type(k_set) :: q_set_uc
        !> Unit-cell (G+q)-grid
        type(Gk_set) :: gq_set_uc
        !> Supercell (G+q)-grid
        type(Gk_set) :: gq_set_sc
        !> Supercell q-grid
        type(k_set) :: q_set_sc
        !> Number of super-cell q-vectors
        integer, parameter :: n_qvecs_sc = 2
        !> Number of super-cell q-vectors
        integer :: n_qvecs_uc
        !> Number of unit-cell (G+q)-vectors per q-vector
        integer, allocatable :: n_gqvecs_uc(:)
        !> Number of unit-cell (G+q)-vectors per q-vector
        integer, allocatable :: n_gqvecs_sc(:)
        !> Supercell dimensions
        integer :: supercell_dims(3)
        !> Running indices
        integer :: i, j
        !> Running indices (G+q)-vectors
        integer :: igq, jgq, iq
        !> Dielectric matrix of the individual component
        !> on the (G+q)-vectors of the heterostructure
        type(dielectric_matrix_type):: eps_sc
        !> Dielectric matrix of the individual component
        !> on the (G+q)-vectors of the heterostructure
        type(dielectric_matrix_type):: eps_sc_ref
        !> Dielectric matrix of the individual component
        !> on the (G+q)-vectors of its sub-unit cell
        type(dielectric_matrix_type):: eps_uc
        complex(dp) :: val
        !> Basis vectors of the supercell reciprocal lattice (columnwise)
        real(dp) :: reciprocal_lattice(3, 3)
        !> Number of spin channels (needed for Gk_set type)
        integer, parameter :: n_spin = 1
        !> Spin channel (needed for Gk_set type)
        integer, parameter :: ispin = 1
        !> Index of Gamma point, i.e. q=0
        integer:: id_gamma_sc
        !> Index of Gamma point, i.e. q=0
        integer  :: id_gamma_uc
        !> List of indices corresponding to zero vector
        integer, allocatable :: indices_zero_q(:)


        ! Setup the basis vectors of the reciprocal supercell lattice
        reciprocal_lattice(:, 1) = (/1._dp, 0._dp, 0._dp/)
        reciprocal_lattice(:, 2) = (/0._dp, 1._dp, 0._dp/)
        reciprocal_lattice(:, 3) = (/0._dp, 0._dp, 2._dp/)

        ! Setup supercell dimensions (2 x 2 x 1)
        supercell_dims = (/2, 2, 1/)

        ! Number of unit cell q-vectors
        n_qvecs_uc = n_qvecs_sc*product(supercell_dims)

        allocate (n_gqvecs_uc(n_qvecs_uc))
        allocate (n_gqvecs_sc(n_qvecs_sc))

        ! Only 1 (G+q)-vector per q-vector for unit cell
        ! (only q-vector itself, i.e. G=0 considered)
        n_gqvecs_uc = 1

        ! Each 4 unit-cell q-vectors fold back to one supercell q-vector
        n_gqvecs_sc = 4

        allocate (q_set_uc%vkc(3, n_qvecs_uc))
        allocate (gq_set_uc%ngk(n_spin, n_qvecs_uc))
        gq_set_uc%ngk(ispin, :) = n_gqvecs_uc

        allocate (gq_set_uc%vgkc(3, maxval(gq_set_uc%ngk), n_spin, n_qvecs_uc))

        q_set_uc%nkpt = n_qvecs_uc

        allocate (q_set_sc%vkc(3, n_qvecs_sc))
        allocate (gq_set_sc%ngk(n_spin, n_qvecs_sc))
        gq_set_sc%ngk(ispin, :) = n_gqvecs_sc

        allocate (gq_set_sc%vgkc(3, maxval(gq_set_sc%ngk), n_spin, n_qvecs_sc))

        q_set_sc%nkpt = n_qvecs_sc

        ! Setup two supercell q-vectors
        q_set_sc%vkc(:, 1) = (/0._dp, 0._dp, 0._dp/)
        q_set_sc%vkc(:, 2) = (/.5_dp, 0._dp, 0._dp/)

        ! Construct 8 unit-cell q-vectors by shifting supercell q-vectors
        ! by reciprocal lattice vectors
        q_set_uc%vkc(:, 1) = (/0._dp, 0._dp, 0._dp/)
        q_set_uc%vkc(:, 2) = (/0.0_dp, 1._dp, 0._dp/)
        q_set_uc%vkc(:, 3) = (/1._dp, 0._dp, 0._dp/)
        q_set_uc%vkc(:, 4) = (/1._dp, 1._dp, 0._dp/)
        q_set_uc%vkc(:, 5) = (/0.5_dp, 0._dp, 0._dp/)
        q_set_uc%vkc(:, 6) = (/.5_dp, 1._dp, 0._dp/)
        q_set_uc%vkc(:, 7) = (/1.5_dp, 0._dp, 0._dp/)
        q_set_uc%vkc(:, 8) = (/1.5_dp, 1._dp, 0._dp/)

        ! Supercell (G+q)-vectors given by unit-cell q-vectors
        do i = 1, gq_set_sc%ngk(ispin, 1)
            gq_set_sc%vgkc(:, i, ispin, 1) = q_set_uc%vkc(:, i)
        end do

        do i = 1, n_gqvecs_sc(2)
            gq_set_sc%vgkc(:, i, ispin, 2) = q_set_uc%vkc(:, n_gqvecs_sc(1) + i)
        end do

        ! Unit-cell (G+q)-vectors given by unit-cell q-vectors (i.e. only G=0)
        do i = 1, n_qvecs_uc
            gq_set_uc%vgkc(:, 1, ispin, i) = q_set_uc%vkc(:, i)
        end do

        ! Find indices of q=0

        if(allocated(indices_zero_q)) deallocate(indices_zero_q)
        indices_zero_q = indices_zero_vectors(q_set_sc%vkc)
        id_gamma_sc = indices_zero_q(1)
        if(allocated(indices_zero_q)) deallocate(indices_zero_q)
        indices_zero_q = indices_zero_vectors(q_set_uc%vkc)
        id_gamma_uc = indices_zero_q(1)



        ! Mock unit-cell dielectric matrix
        call eps_uc%init(gq_set_uc%ngk(ispin, :), id_gamma_uc)

        do iq = 1, q_set_uc%nkpt
            do igq = 1, gq_set_uc%ngk(ispin, iq)
                val = cmplx(norm2(gq_set_uc%vgkc(:, igq, ispin, iq)), 0._dp, dp)
                eps_uc%body(igq, igq, iq) = val
            end do
        end do

        eps_uc%head = zone
        eps_uc%wings = zone

        call eps_sc%init(gq_set_sc%ngk(ispin, :), id_gamma_sc)

        ! Expand unit-cell dielectric matrix to supercell
        call expand_epsilon(gq_set_uc, q_set_uc, gq_set_sc, q_set_sc, reciprocal_lattice, &
                            supercell_dims, eps_uc, eps_sc)

        ! Create expected expanded supercell matrix
        call eps_sc_ref%init(gq_set_sc%ngk(ispin, :), id_gamma_sc)

        do iq = 1, q_set_sc%nkpt
            do igq = 1, gq_set_sc%ngk(ispin, iq)
                val = cmplx(norm2(gq_set_sc%vgkc(:, igq, ispin, iq)), 0._dp, dp)
                eps_sc_ref%body(igq, igq, iq) = val
            end do
        end do

        call test_report%assert(all_close(eps_sc_ref%body, eps_sc%body), &
                                'Test test_expand_epsilon: &
                                Expected: Last 4 unit-cell q-vectors fold &
                                back to second supercell q-vector.')
        call eps_sc%delete
        call eps_sc_ref%delete
        call eps_uc%delete

    end subroutine

    !> Test fold_uc_to_sc_qgrid,
    !> i.e. finding the indices of all unit-cell q-vectors
    !> folding back to a supercell q-vector. For testing, a supercell
    !> and a unit-cell q-grid are set up, where the unit-cell q-vectors folding
    !> back to a given supercell q-vector are constructed as the sum of the
    !> supercell q-vector and supercell G-vectors.
    !> Note that the term 'supercell' refers to real space, such that
    !> the supercell Brillouin zone is smaller than the unit cell BZ.
    !> For more information about the connection between the two q-grids see
    !> the routine documentation.
    subroutine test_fold_uc_to_sc_qgrid(test_report)
        use math_utils, only: all_close

        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report
        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Unit-cell q-vectors
        real(dp), allocatable :: q_vecs_uc(:, :)
        !> Super-cell q-vectors
        real(dp), allocatable :: q_vecs_sc(:, :)
        !> Supercell dimensions
        integer :: supercell_dims(3)
        !> Number of super-cell q-vectors
        integer, parameter :: n_qvecs_sc = 2
        !> Number of super-cell q-vectors
        integer :: n_qvecs_uc
        !> Basis vectors of the supercell reciprocal lattice (columnwise)
        real(dp) :: reciprocal_lattice(3, 3)

        !> Indices of all unit-cell q-vectors folding back to uc-vector
        integer, allocatable :: map_sc2uc(:)

        supercell_dims = (/2, 2, 1/)

        n_qvecs_uc = n_qvecs_sc*product(supercell_dims)

        allocate (q_vecs_uc(3, n_qvecs_uc))
        allocate (q_vecs_sc(3, n_qvecs_sc))

        ! Setup two supercell q-vectors
        q_vecs_sc(:, 1) = (/0._dp, 0._dp, 0._dp/)
        q_vecs_sc(:, 2) = (/.5_dp, 0._dp, 0._dp/)

        ! Setup the basis vectors of the reciprocal supercell lattice
        reciprocal_lattice(:, 1) = (/1._dp, 0._dp, 0._dp/)
        reciprocal_lattice(:, 2) = (/0._dp, 1._dp, 0._dp/)
        reciprocal_lattice(:, 3) = (/0._dp, 0._dp, 2._dp/)

        ! Setup 4 unit-cell q-vectors folding back to first supercell q-vector
        q_vecs_uc(:, 1) = (/0._dp, 0._dp, 0._dp/)
        q_vecs_uc(:, 2) = (/0.0_dp, 1._dp, 0._dp/)
        q_vecs_uc(:, 3) = (/1._dp, 0._dp, 0._dp/)
        q_vecs_uc(:, 4) = (/1._dp, 1._dp, 0._dp/)

        ! Setup 4 unit-cell q-vectors folding back to second supercell q-vector
        q_vecs_uc(:, 5) = (/0.5_dp, 0._dp, 0._dp/)
        q_vecs_uc(:, 6) = (/.5_dp, 1._dp, 0._dp/)
        q_vecs_uc(:, 7) = (/1.5_dp, 0._dp, 0._dp/)
        q_vecs_uc(:, 8) = (/1.5_dp, 1._dp, 0._dp/)

        ! Find unit cell q-vectors folding to first supercell vector
        map_sc2uc = fold_uc_to_sc_qgrid(supercell_dims, reciprocal_lattice, &
                                        q_vecs_uc, q_vecs_sc(:, 1))

        call test_report%assert(all(map_sc2uc == (/1, 2, 3, 4/)), &
                                'Test sc2uc_map: &
                                Expected: First 4 unit-cell q-vectors fold &
                                back to first supercell q-vector.')

        deallocate (map_sc2uc)

        map_sc2uc = fold_uc_to_sc_qgrid(supercell_dims, reciprocal_lattice, &
                                        q_vecs_uc, q_vecs_sc(:, 2))

        call test_report%assert(all(map_sc2uc == (/5, 6, 7, 8/)), &
                                'Test sc2uc_map: &
                                Expected: Last 4 unit-cell q-vectors fold &
                                back to second supercell q-vector.')

    end subroutine

end module expand_add_eps_tests
