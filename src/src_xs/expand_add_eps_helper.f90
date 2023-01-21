!> Routines needed in the addition and expansion approach
!> implemented in expand_add_eps.f90
module expand_add_eps_helper

    use modmpi, only: mpiglobal, terminate_mpi_env
    use precision, only: sp, dp
    

    implicit none

    private

    public ::   &
              setup_gq_vecs_from_file, setup_q_vecs_from_file, &
              setup_gq_vecs_from_data, setup_q_vecs_from_data, &
               fold_uc_to_sc_qgrid

    real(dp), parameter  :: tol_default = 1e-9_dp


contains

    !> Reads q-vectors from file and stores them in the respective
    !> member of the 'Gk_set' type.
    subroutine setup_q_vecs_from_file(filename, n_qvecs, q_set)
        use mod_kpointset, only: k_set
        use xs_file_interface, only: read_cartesian_coordinates
        use modmpi, only: mpiglobal
        use exciting_mpi, only: xmpi_bcast

        character(*), intent(in) :: filename
        class(k_set), intent(out) :: q_set
        !> Number of q-vectors
        integer, intent(in) :: n_qvecs

        allocate (q_set%vkc(3, n_qvecs))

        if (mpiglobal%is_root) then
            call read_cartesian_coordinates(trim(filename), &
                                            q_set%vkc)
        end if
        call xmpi_bcast(mpiglobal, q_set%vkc)

        q_set%nkpt = n_qvecs

    end subroutine

    !> Stores given data in the respective member of the 'k_set' type.
    subroutine setup_q_vecs_from_data(q_vecs, q_set)
        use mod_kpointset, only: k_set, delete_k_vectors

        !> q-vectors in cartesian coordinates
        real(dp), intent(in) ::  q_vecs(:, :)
        class(k_set), intent(out) :: q_set
        !> Number of q-vectors
        integer :: n_qvecs

        n_qvecs = size(q_vecs, dim=2)
        allocate (q_set%vkc(3, n_qvecs), source=q_vecs)
        q_set%nkpt = n_qvecs

    end subroutine

    !> Reads (G+q)-vectors from file and stores them in the respective
    !> member of the 'Gk_set' type.
    subroutine setup_gq_vecs_from_file(gqvecs_dirname, qvecs_fname, n_spin, n_qvecs, gq_set)
        use mod_kpointset, only: Gk_set
        use xs_file_interface, only: read_n_gq, read_gq_vectors
        use modmpi, only: mpiglobal
        use exciting_mpi, only: xmpi_bcast

        !> Name of file containing q-vectors
        character(*), intent(in) :: qvecs_fname
        !> Name of directory containing (G+q)-vectors
        character(*), intent(in) :: gqvecs_dirname
        class(Gk_set), intent(out) :: gq_set
        !> Number of q-vectors
        integer, intent(in) :: n_qvecs
        !> Number of spin channels (needed for Gk_set type)
        integer, intent(in) :: n_spin
        !> Spin channel (needed for Gk_set type)
        integer, parameter :: ispin = 1

        allocate (gq_set%ngk(n_spin, n_qvecs))

        if (mpiglobal%is_root) then
            call read_n_gq(trim(qvecs_fname), gq_set%ngk(ispin, :))
        end if
        call xmpi_bcast(mpiglobal, gq_set%ngk(ispin, :))
        allocate (gq_set%vgkc(3, maxval(gq_set%ngk), n_spin, n_qvecs))

        if (mpiglobal%is_root) then
            call read_gq_vectors(trim(gqvecs_dirname), &
                                 gq_set%ngk(ispin, :), gq_set%vgkc(:, :, ispin, :))
        end if

        call xmpi_bcast(mpiglobal, gq_set%vgkc(:, :, ispin, :))

        gq_set%ngkmax = maxval(gq_set%ngk)

    end subroutine

    !> Stores given data in the respective member of the 'Gk_set' type.
    subroutine setup_gq_vecs_from_data(gq_vecs, n_gqvecs, n_spin, gq_set)
        use mod_kpointset, only: Gk_set
        use xs_file_interface, only: read_n_gq, read_gq_vectors
        use modmpi, only: mpiglobal
        use exciting_mpi, only: xmpi_bcast
        !> (G+q)-vectors in cartesian coordinates
        real(dp), intent(in) :: gq_vecs(:, :, :)
        !> Number of (G+q)-vectors per q-vector
        integer, intent(in) :: n_gqvecs(:)
        type(Gk_set), intent(out) :: gq_set

        !> Number of spin channels (needed for Gk_set type)
        integer, intent(in) :: n_spin
        !> Spin channel (needed for Gk_set type)
        integer, parameter :: ispin = 1
        !> Maximum number of (G+q)-vectors
        integer :: ngq_max
        !> Number of q-vectors
        integer :: n_qvecs

        n_qvecs = size(n_gqvecs)
        ngq_max = maxval(n_gqvecs)

        allocate (gq_set%ngk(n_spin, n_qvecs))
        allocate (gq_set%vgkc(3, ngq_max, n_spin, n_qvecs))

        gq_set%vgkc(:, :, ispin, :) = gq_vecs
        gq_set%ngk(ispin, :) = n_gqvecs
        gq_set%ngkmax = ngq_max

    end subroutine

    


    !> Creates a map for a given q-vector from a m x n x 1 supercell
    !> which contains the indices of all unit-cell q-vectors folding
    !> to the supercell q-vector.
    !> Note that the terms 'supercell' and unit cell'
    !> refer to real space, such that
    !> the supercell Brillouin zone is smaller than the unit cell BZ.
    !> If the q-density is equivalent, multiple q-vectors from the UC
    !> fold back to the SC. Folding means that a
    !> unit-cell q-vector can be expressed as the sum of the supercell
    !> q-vector and the basis vectors of the supercell reciprocal lattice, i.e.
    !>
    !> \[
    !>    \mathbf{q}_{uc} = \mathbf{q}_{sc} + i * \mathbf{b}^{(1)}_{sc}
    !>                                      + j * \mathbf{b}^{(2)}_{sc},
    !> \]
    !>
    !> for some \( i, j \) with
    !>
    !>    \[
    !>      i \in \{0, 1, ..., m - 1\}, j \in \{0, 1, ..., n - 1\}.
    !>                                                              \]
    !>
    !> The number of unit-cell q-vectors folding to one supercell q-vector
    !> is determined by the dimension of the supercell: m * n
    !> unit-cell vectors fold back to one supercell q-vector.
    function fold_uc_to_sc_qgrid(supercell_dims, reciprocal_lattice, qvecs_uc, qvec_sc) result(map)
        use math_utils, only: all_close
        use grid_utils, only: index_column_vector_in_array

        !> All unit-cell q-vectors
        real(dp), intent(in) :: qvecs_uc(:, :)
        !> One supercell q-vector
        real(dp), intent(in) :: qvec_sc(:)
        !> Supercell reciprocal lattice vectors, stored columnwise
        real(dp), intent(in) :: reciprocal_lattice(3, 3)
        !> Supercell dimensions (m x n x 1)
        integer, intent(in) :: supercell_dims(3)
        !> Indices of all unit-cell q-vectors folding
        !> to supercell q-vector
        integer, allocatable :: map(:)
        !> Running index unit-cell q-vector
        integer :: iq_uc
        !> Index counting the folding q-vectors
        integer :: i_fold
        !> Running indices over supercell dimensions
        integer ::  i, j
        !> Number of unit-cell q-vectors
        integer :: n_qvecs_uc
        !> G-vector from the supercell
        real(dp) :: g_vec_sc(3)
        !> Unit-cell q-vector obtained by adding a
        !> G-vector from the supercell to a supercell q-vector
        real(dp) :: q_vec(3)

        n_qvecs_uc = size(qvecs_uc, dim=2)
        allocate (map(product(supercell_dims)))

        ! Third component is equivalent for uc and sc
        q_vec(3) = reciprocal_lattice(3, 3)
        g_vec_sc(3) = 0._dp

        ! Construct all unit cell q-vectors folding to supercell q-vector
        ! and find their indices.
        i_fold = 0
        do i = 0, supercell_dims(1) - 1
            do j = 0, supercell_dims(2) - 1
                g_vec_sc(1:2) = matmul(reciprocal_lattice(1:2, 1:2), (/i, j/))

                q_vec = qvec_sc + g_vec_sc

                ! Find index in set of unit-cell q-vectors
                iq_uc = index_column_vector_in_array(q_vec, qvecs_uc, tol=1e-6_dp)

                i_fold = i_fold + 1
                map(i_fold) = iq_uc
            end do
        end do

    end function

end module
