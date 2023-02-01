!> Module for the "XAF"-method (without F="full wavefunction")
!> presented first in J. Chem. Theory Comput. 2019, 15, 3824−383.
!> A similar procedure is presented in
!> J. Chem. Theory Comput. 2019, 15, 4218−4227,
!> however with a much clearer explanation.
!> The objective is to compute the dielectric matrix of
!> a bilayer heterostructure which can be partitioned into individual systems
!> without the creation of dangling bonds.
!> The two main steps are:
!>
!> "X" = Expansion:
!> The smallest sub-unit cell for each separate component is
!> identified, and the dielectric matrix of the
!> component is computed using this sub-unit cell.
!> The dielectric matrix of the same component in the supercell used in the
!> heterostructure is obtained by an expansion (unfolding).
!>
!> The calculation of the individual dielectric matrices
!> has to be done in advance, i.e. before calling this routine,
!> and the dielectric matrices of the two systems have to be present in the
!> run directory in the sub-directories EPS0_1 and EPS0_2, respectively.
!> Additionally, the (G+q)-vectors of the sub-unit cells used in the
!> calculation of the individual components have to be present in the
!> directories GQPOINTS_1 and GQPOINTS_2, respectively.
!>
!> "A" = Add:
!> The DM of the of the heterostructure is obtained by adding
!> the polarizabilities (which can be obtained through the DM)
!> of the two systems.
module expand_add_eps
    use precision, only: sp, dp
    use asserts, only: assert
    use expand_add_eps_helper, only: setup_gq_vecs_from_file, &
    setup_q_vecs_from_file, &
    setup_gq_vecs_from_data, &
    setup_q_vecs_from_data, fold_uc_to_sc_qgrid


    implicit none

    private
    public :: expand_add_eps_launcher, expand_epsilon

contains

    !> Launches the XAF method: Initializes all global exciting variables
    !> of the heterostrucuture needed for
    !> the calculation such as the reciprocal lattice vectors
    !> (bvec), the (G+q)-vectors in cartesian/lattice coordinates (vgqc/vgql),
    !> the q-vectors in cartesian/lattice coordinates (vqc/vql), and the number
    !> of (G+q)-vectors per q-vector (ngq), and passes them to the main routine.
    subroutine expand_add_eps_launcher
        use modmain, only: bvec
        use mod_qpoint, only: vqc, vql, nqpt
        use modxs, only: vgqc, ngq, vgql
        use modinput, only: input
        use modmpi, only: mpiglobal
        use os_utils, only: make_directory_command 
        use xs_file_interface, only: write_gq_vectors, write_q_vectors
        !> System command
        character(256) :: syscommand
        !> Number of components (fixed to 2)
        integer, parameter :: n_components = 2
        !> Container for supercell dimensions of the two individual components
        integer :: supercell_dims(3, n_components)
        !> Name of directory containing (G+q)-vectors for the heterostructure
        character(256), parameter :: gqvecs_dirname = 'GQPOINTS'
        !> Name of directory containing q-vectors for the heterostructure
        character(256), parameter :: qvecs_fname = 'QPOINTS_SCR.OUT'
        !> Name of directory containing (G+q)-vectors for the heterostructure
        character(256), parameter :: eps0_dirname = 'EPS0'
        !> Running index q-vectors
        integer :: iq
        !> OS command
        character(256) :: os_command

        ! Initialise (reciprocal) lattice parameters and (G+q)-grids
        call init0
        call init1
        call init2

        ! Making folder for heterostructure DMs and (G+q)-vectors
        if (mpiglobal%is_root) then

            os_command = make_directory_command(eps0_dirname)
            call system(trim(adjustl(os_command)))
            os_command = make_directory_command(gqvecs_dirname)
            call system(trim(adjustl(os_command)))


            ! Write grid of heterostructure to file
            call write_q_vectors(qvecs_fname, vql, vqc, ngq)
            do iq = 1, nqpt
                call write_gq_vectors(gqvecs_dirname, iq, vgql(:, :, iq), &
                                      vgqc(:, :, iq), ngq(iq))
            end do
        end if
        supercell_dims(:, 1) = input%xs%expand_eps%supercell_1
        supercell_dims(:, 2) = input%xs%expand_eps%supercell_2

        call expand_add_epsilon_main(bvec, supercell_dims, &
                                     vgqc, vqc, ngq)

    end subroutine

    !> Main routine for expansion and addition of the
    !> dielectric matrix (DM). Responsible for all I/O, reads:
    !>       i) (G+q)-vectors and q-vectors of the sub-unit cells
    !>            of the individual components,
    !>       ii) Dielectric matrix of the sub-unit cells
    !>            of the individual components.
    !> Data is passed to expand_epsilon, which generates the DM of the
    !> individual components on the (G+q)-grid of the heterostructure.
    !> Finally, the DM of the full heterostructure can be obtained
    !> by adding the polarizabilites of the two systems.
    !> For \( \mathbf{q} = 0 \), the elements
    !> \(  \varepsilon_{00}(0) \) (dimension 3 x 3),
    !> \(  \varepsilon_{G0}(0) \) and \(  \varepsilon_{0G'}(0) \)
    !> (dimension n_gqvecs x 2 x 3) are referred to as head and wings,
    !> respectively.
    subroutine expand_add_epsilon_main(bvec, &
                                       supercell_dims, gq_vecs_hs, &
                                       q_vecs_hs, &
                                       n_gq_vecs_hs)
        use diel_mat_type, only: dielectric_matrix_type, compare_dielectric_matrix, write_dielectric_matrix
 
        use xs_file_interface, only: read_gq_vectors, &
                                     read_cartesian_coordinates, read_n_gq
        use constants, only: zzero
        use mod_kpointset, only: Gk_set, k_set, &
                                 delete_Gk_vectors, delete_k_vectors
        use grid_utils, only: indices_zero_vectors

        !> Reciprocal lattice vectors, stored column-wise
        real(dp), intent(in) :: bvec(3, 3)
        !> (G+q)-vecs of the heterostructure
        real(dp), intent(in) :: gq_vecs_hs(:, :, :)
        !> (G+q)-vecs of the heterostructure
        real(dp), intent(in) :: q_vecs_hs(:, :)
        !> (G+q)-vecs of the heterostructure
        integer, intent(in) :: n_gq_vecs_hs(:)

        !> Index of first individual component for container types
        integer, parameter :: ind_1 = 1
        !> Index of first individual component for container types
        integer, parameter :: ind_2 = 2
        !> Index of heterostructure for container types
        integer, parameter :: hs = 3
        !> Supercell dimensions of first component
        integer, intent(in) :: supercell_dims(3, 2)

        !> (G+q)-grids of the three systems
        type(Gk_set) :: gq_grids(3)

        !> DM of the two individual systems on their respective (G+q)-grids
        type(dielectric_matrix_type):: eps_on_ind_grid(2)
        !> q-grids for the three systems (two ind. and heterostructure)
        type(k_set) :: q_grids(3)

        !> DM of the three systems on (G+q)-grid of the heterostructure
        type(dielectric_matrix_type):: eps_on_hs_grid(3)

        !> Number of q-vectors for the three systems
        integer :: n_qvecs(3)

        !> Name of files containing q-vectors of the three systems
        character(256), parameter :: qvecs_fnames(3) = (/character(len=256) :: &
                                                         'QPOINTS_SCR_1.OUT', &
                                                         'QPOINTS_SCR_2.OUT', &
                                                         'QPOINTS_SCR.OUT'/)
        !> Name of directories containing (G+q)-vectors for the three systems
        character(256), parameter :: gqvecs_dirnames(3) = (/character(len=256):: &
                                                            'GQPOINTS_1', &
                                                            'GQPOINTS_2', &
                                                            'GQPOINTS'/)
        !> Name of directories containing diel. matrices for the three systems
        character(256), parameter :: eps0_dirnames(3) = (/character(len=256) :: &
                                                          'EPS0_1', &
                                                          'EPS0_2', &
                                                          'EPS0'/)

        !> Index of Gamma point (q=0) for the q-grid of the heterostructure
        integer :: id_gamma

        !> Dummy index
        integer :: i
        character(256) :: strnum, fname
        !> Number of spin channels (needed for Gk_set type)
        integer, parameter :: n_spin = 1
        !> Spin channel (needed for Gk_set type)
        integer, parameter :: ispin = 1
        !> Running index systems
        integer :: i_system
        !> Index of q = 0
        integer, allocatable :: indices_zero_q(:)

        n_qvecs(hs) = size(n_gq_vecs_hs)

        ! Read (G+q)-vectors of the two individual systems
        do i_system = 1, 2
            n_qvecs(i_system) = n_qvecs(hs)*product(supercell_dims(:, i_system))
            call setup_q_vecs_from_file(qvecs_fnames(i_system), &
                                        n_qvecs(i_system), q_grids(i_system))

            call setup_gq_vecs_from_file(gqvecs_dirnames(i_system), &
                                         qvecs_fnames(i_system), &
                                         n_spin, n_qvecs(i_system), &
                                         gq_grids(i_system))

        end do

        ! Setup (G+q)-vectors for HS from input data (modxs variables)
        call setup_q_vecs_from_data(q_vecs_hs, q_grids(3))

        call setup_gq_vecs_from_data(gq_vecs_hs, &
                                     n_gq_vecs_hs, n_spin, &
                                     gq_grids(3))

        ! Index of q=0 in q-grid of the heterostructure
        indices_zero_q = indices_zero_vectors(q_grids(hs)%vkc)
        id_gamma = indices_zero_q(1)

        ! Allocate dielectric matrices on (G+q)-grid of the heterostructure
        do i_system = 1, 3
            call eps_on_hs_grid(i_system)%init(gq_grids(hs)%ngk(ispin, :), &
                                               id_gamma)
        end do

        ! Expand DMs of the ind. components to (G+q)-grid of the heterostructure
        do i_system = 1, 2

            call eps_on_ind_grid(i_system)%init(gq_grids(i_system)%ngk(ispin, :), id_gamma)

            ! Initialize and read sub-unit cell DM of the individual component
            call eps_on_ind_grid(i_system)%setup_from_file(eps0_dirnames(i_system), &
                                  & gq_grids(i_system)%ngk(ispin, :), &
                                    q_grids(i_system)%vkc)

            ! Obtain DM of ind. component on (G+q)-grid of the heterostructure
            call expand_epsilon(gq_grids(i_system), q_grids(i_system), &
                                gq_grids(hs), q_grids(hs), &
                                bvec, supercell_dims(:, i_system), eps_on_ind_grid(i_system), eps_on_hs_grid(i_system))

            call eps_on_ind_grid(i_system)%delete
        end do

        ! Obtain DM of heterostructure by adding individual polarizabilities
        eps_on_hs_grid(hs) = add_dielectric_matrices(eps_on_hs_grid(ind_1), &
                                                  eps_on_hs_grid(ind_2))

        ! Write dielectric matrix of heterostructure to file
        call  eps_on_hs_grid(hs)%write_to_file(eps0_dirnames(hs), &
                                     gq_grids(hs)%ngk(ispin, :), &
                                     q_grids(hs)%vkc)

        ! Free memory
        do i_system = 1, 3
            call eps_on_hs_grid(i_system)%delete
            call delete_Gk_vectors(gq_grids(i_system))
            call delete_k_vectors(q_grids(i_system))
        end do
        
    end subroutine

    !> Given the two dielectric matrices (DMs),
    !> \( \varepsilon^{1,2}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) \),
    !> this function computes the dielectric matrix which is obtained when
    !> the polarizabilites \( \chi^{1,2}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) \)
    !> of the two systems are added.
    !> As polarizability and DM are connected via
    !>
    !> \[
    !>   \varepsilon_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) =
    !>       \delta_{\mathbf{G}\mathbf{G}'}
    !>     - v_\mathbf{G}(\mathbf{q})\chi_{\mathbf{G}\mathbf{G}'}(\mathbf{q}).
    !> \],
    !>
    !> the output dielectric matrix is obtained as
    !>
    !> \[
    !>    \varepsilon_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) =
    !>        \varepsilon^1_{\mathbf{G}\mathbf{G}'}(\mathbf{q})
    !>      + \varepsilon^2_{\mathbf{G}\mathbf{G}'}(\mathbf{q})
    !>      - \delta_{\mathbf{G}\mathbf{G}'}
    !> \].
    function add_dielectric_matrices(eps_1, eps_2) result(eps_out)
        use diel_mat_type, only: dielectric_matrix_type        
        use math_utils, only: identity_complex_dp

        type(dielectric_matrix_type):: eps_1, eps_2, eps_out

        integer :: n_qvecs, n_gqmax, iq

        n_qvecs = size(eps_1%body, dim=3)
        n_gqmax = size(eps_1%body, dim=2)
        ! Add contribution of individual component
        eps_out%head = eps_1%head + eps_2%head
        eps_out%wings = eps_1%wings + eps_2%wings
        eps_out%body = eps_1%body + eps_2%body

        ! Subtract 1 from the diagonals of head and body
        ! (polarizability is additive and not the DM, see routine documentation)
        eps_out%head = eps_out%head - identity_complex_dp(3)

        do iq = 1, n_qvecs

            eps_out%body(:, :, iq) = eps_out%body(:, :, iq) &
                                     - identity_complex_dp(n_gqmax)
        end do

    end function

    !> Expands the dielectric matrix (DM) computed for a unit cell,
    !>
    !> \[
    !>      \varepsilon_{\mathbf{G_{uc}G_{uc}'}}(\mathbf{q}_{uc}),
    !> \]
    !>
    !> to the DM corresponding to a supercell,
    !>
    !> \[
    !>      \varepsilon_{\mathbf{G}_{sc}\mathbf{G}_{sc}'}(\mathbf{q}_{sc}).
    !> \]
    !>
    !> The expansion is done by identifying the elements of the unit-cell
    !> dielectric matrix with the elements of the supercell, i.e. each element
    !> of the unit-cell gets mapped to one of the supercell:
    !>
    !>\[
    !>   \varepsilon_{\mathbf{G}_{uc}\mathbf{G}_{uc}'}(\mathbf{q}_{uc})
    !>                  \rightarrow
    !>   \varepsilon_{\mathbf{G}_{sc}\mathbf{G}_{sc}'}(\mathbf{q}_{sc})
    !>                                                                  \]
    !> with
    !>
    !> \[
    !>   \mathbf{q}_{uc} + \mathbf{G}_{uc} = \mathbf{q}_{sc} + \mathbf{G}_{sc},
    !>   \mathbf{q}_{uc} + \mathbf{G}'_{uc} = \mathbf{q}_{sc} + \mathbf{G}'_{sc}
    !>                                                                        \]
    !>
    !> Not all elements of the supercell dielectric matrix will be targeted by
    !> this mapping procedure, all these elements have to be zero such that
    !> the amount of information in the unit cell matrix is equivalent
    !> to the one in the supercell matrix. Using this mapping procedure,
    !> the dielectric matrix of a given supercell q-vector contains
    !> contributions from different unit-cell q-vectors. These vectors
    !> are referred to as the ones folding to  the supercell q-vector
    !>  and can be determined  on the fly during the expansion.
    subroutine expand_epsilon(gq_set_uc, q_set_uc, gq_set_sc, q_set_sc, &
                              reciprocal_lattice, &
                              supercell, eps_uc, eps_sc)

        use grid_utils, only: index_column_vector_in_array
        use constants, only: zzero
        use math_utils, only: identity_complex_dp
        use modmpi, only: mpiglobal, distribute_loop, &
                          mpi_allgatherv_ifc, terminate_if_false
        use mod_kpointset, only: Gk_set, k_set
        use diel_mat_type, only: dielectric_matrix_type
        use grid_utils, only: indices_zero_vectors
        use exciting_mpi, only: xmpi_bcast

        !> (G+q)-vectors of the unit cell in cartesian coordinates
        type(Gk_set), intent(inout)  :: gq_set_uc
        !> q-vectors of the unit cell in cartesian coordinates
        type(k_set), intent(inout)  :: q_set_uc
        !> (G+q)-vectors of the supercell in cartesian coordinates
        type(Gk_set), intent(inout)  :: gq_set_sc
        !> q-vectors of the supercell in cartesian coordinates
        type(k_set), intent(inout)  :: q_set_sc

        !> Reciprocal lattice vectors, stored columnwise
        real(dp), intent(in) :: reciprocal_lattice(3, 3)
        !> Dimension of the supercell
        integer, intent(in) ::  supercell(:)

        !> Unit-cell dielectric matrix
        type(dielectric_matrix_type), intent(in)  :: eps_uc
        !> Expanded supercell dielectric matrix
        type(dielectric_matrix_type), intent(inout)  :: eps_sc

        !> Index of (G+q)-vectors in unit cell
        integer :: igq_uc
        !> Index of (G+q)-vectors in unit cell
        integer :: jgq_uc
        !> Index of q-vectors in unit cell
        integer :: iq_uc
        !> Index of (G+q)-vector in supercell
        integer :: igq_sc
        !> Index of q-vector in supercell
        integer :: iq_sc
        !> Index of q-vector in supercell
        integer :: jq_sc
        !> Index of (G+q)-vector in supercell
        integer :: jgq_sc
        !> Dummy index
        integer :: i
        !> Indices of all unit-cell q-vectors folding back to supercell q-vector
        integer, allocatable :: map_uc2sc(:)
        !> Spin channel (needed for Gk_set type)
        integer, parameter :: ispin = 1
        !> Index of Gamma point, i.e. q=0
        integer, allocatable :: id_gamma_uc
        !> First q-vector for current rank for mpi distribution
        integer :: iq_start
        !> First q-vector for current rank for mpi distribution
        integer :: iq_end
        !> Indices of all zero vectors
        integer, allocatable ::  indices_zero_q(:)


        ! All elements not present in the unit cell will not be touched
        eps_sc%wings = zzero
        eps_sc%body = zzero

        indices_zero_q = indices_zero_vectors(q_set_uc%vkc)
        id_gamma_uc = indices_zero_q(1)

        if (mpiglobal%is_root) then

            ! Super-cell head equals unit-cell head
            eps_sc%head = eps_uc%head
            
            ! Super-cell wings can only have contributions from unit-cell wings
            do igq_uc = 1, gq_set_uc%ngk(ispin, id_gamma_uc)
                igq_sc = index_column_vector_in_array( &
                         gq_set_uc%vgkc(:, igq_uc, ispin, id_gamma_uc), &
                         gq_set_sc%vgkc(:, :, ispin, 1), tol=1e-6_dp)
                call terminate_if_false(igq_sc > 0, &
                    & message='Error expand_epsilon:&
                    &Supercell vector not in unit-cell grid.')

                eps_sc%wings(igq_sc, :, :) = eps_uc%wings(igq_uc, :, :)
            end do

        end if
        call xmpi_bcast(mpiglobal, eps_sc%head)
        call xmpi_bcast(mpiglobal, eps_sc%wings)

        call distribute_loop(mpiglobal,q_set_sc%nkpt,iq_start,iq_end )
        do iq_sc = iq_start, iq_end

            ! Find indices of all unit-cell q-vectors folding
            ! back to supercell q-vector with index iq_sc
            map_uc2sc = fold_uc_to_sc_qgrid( &
                        supercell, reciprocal_lattice, &
                        q_set_uc%vkc, q_set_sc%vkc(:, iq_sc))

            do i = 1, size(map_uc2sc)
                iq_uc = map_uc2sc(i)
                do igq_uc = 1, gq_set_uc%ngk(ispin, iq_uc)
                    do jgq_uc = 1, gq_set_uc%ngk(ispin, iq_uc)

                        ! Find index of super-cell (G+q)-vector
                        ! corresponding to (G+q)-vector from the unit cell
                        igq_sc = index_column_vector_in_array( &
                                 gq_set_uc%vgkc(:, igq_uc, ispin, iq_uc), &
                                 gq_set_sc%vgkc(:, :, ispin, iq_sc), &
                                 tol=1e-6_dp)

                        ! Find index of super-cell (G'+q)-vector
                        ! corresponding to (G'+q)-vector from the unit cell
                        jgq_sc = index_column_vector_in_array( &
                                &   gq_set_uc%vgkc(:, jgq_uc, ispin, iq_uc), &
                                    gq_set_sc%vgkc(:, :, ispin, iq_sc), &
                                    tol=1e-6_dp)

                        call terminate_if_false(igq_sc > 0 .and. jgq_sc > 0, &
                            & message='expand_epsilon:&
                            &Supercell vector not in unit-cell grid.')

                        eps_sc%body(igq_sc, jgq_sc, iq_sc) &
                            = eps_uc%body(igq_uc, jgq_uc, iq_uc)
                    end do
                end do
            end do
        end do

        ! Gather body on all processes
        call mpi_allgatherv_ifc(set=q_set_sc%nkpt, &
                                rlen=size(eps_sc%body, dim=1)**2, &
                                zbuf=eps_sc%body, inplace=.true., &
                                comm=mpiglobal)

    end subroutine

end module
