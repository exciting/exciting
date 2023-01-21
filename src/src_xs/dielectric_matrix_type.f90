module diel_mat_type
    use precision, only:  dp

    implicit none

    real(dp), parameter  :: tol_default = 1e-9_dp
    
    !> Dielectric matrix type
    Type dielectric_matrix_type
        !> Head of the dielectric matrix (only present for q=0)
        complex(dp), allocatable  :: head(:, :)
        !> Wings of the dielectric matrix (only present for q=0)
        complex(dp), allocatable  :: wings(:, :, :)
        !> Body of the dielectric matrix for all q-vectors
        complex(dp), allocatable  :: body(:, :, :)

    contains
        !> Initialise by parsing arguments from the command line
        procedure :: init => allocate_dielectric_matrix
        procedure :: delete => deallocate_dielectric_matrix
        procedure :: setup_from_file => setup_dielectric_matrix_from_file
        procedure :: write_to_file => write_dielectric_matrix


    End Type dielectric_matrix_type
contains

    !> Allocates all members of the dielectric_matrix_type
    subroutine allocate_dielectric_matrix(this, n_gqvecs, id_gamma)

        use constants, only: zzero

        class(dielectric_matrix_type), intent(inout) :: this
        !> Number of (G+q)-vectors per q-vector
        integer, intent(in) :: n_gqvecs(:)
        !> Index of q=0
        integer, intent(in) :: id_gamma
        !> Maximum number of (G+q)-vectors
        integer :: n_gqmax
        !> Number of (G+q)-vectors for q=0
        integer:: n_gq_gamma
        
        !> Number of q-vectors
        integer :: n_qvecs

        n_qvecs = size(n_gqvecs)
        n_gq_gamma = n_gqvecs(id_gamma)
        n_gqmax = maxval(n_gqvecs)

        allocate (this%head(3, 3))
        allocate (this%body(n_gqmax, n_gqmax, n_qvecs), source=zzero)
        allocate (this%wings(n_gq_gamma, 2, 3), source=zzero)

    end subroutine

    !> Deallocates all members of the dielectric_matrix_type
    subroutine deallocate_dielectric_matrix(this)

        use constants, only: zzero

        class(dielectric_matrix_type), intent(inout) :: this

        if (allocated(this%head)) deallocate (this%head)
        if (allocated(this%wings)) deallocate (this%wings)
        if (allocated(this%body)) deallocate (this%body)

    end subroutine

    !> A wrapper for geteps0_finite_q and  geteps0_zero_q:
    !> Reads the dielectric matrix for all q-vectors employed in a BSE
    !> calculation from binary files and stores it in the corresponding
    !> arrays of the dielectric matrix type. Files read for zero frequency.
    !> The directory name has to be specified, while the filenames
    !> are expected to be EPS0_QXXXXX.OUT, where XXXXX is the q-point index
    !> with leading zeros.
    subroutine setup_dielectric_matrix_from_file(this, eps0_dirname, n_gqvecs, &
                                                 qvecs)
        use modmpi, only: mpiglobal, mpi_allgatherv_ifc, terminate_if_false
        use exciting_mpi, only: xmpi_bcast
        use math_utils, only: all_zero
        use constants, only: zzero
        use putgeteps0, only: geteps0_finite_q, geteps0_zero_q
        use grid_utils, only: indices_zero_vectors, indices_finite_vectors

        !> Directory where dielectric matrices are stored
        character(*), intent(in) :: eps0_dirname
        !> Number of (G+q)-vectors for all q-vectors
        integer, intent(in) :: n_gqvecs(:)
        !> q-vectors
        real(dp), intent(in) :: qvecs(:, :)
        !> Dielectric matrix
        class(dielectric_matrix_type), intent(inout) :: this

        !> Number of q-vectors
        integer :: n_qvecs
        !> Maximum number of (G+q)-vectors
        integer :: ngq_max
        !> Running index (G+q)-vectors in unit cell
        integer :: iq, i
        !> Name of a file
        character(256) :: fname
        !> string for filename
        character(256) :: strnum
        !> Indices of finite q
        integer, allocatable :: indices_finite_q(:)
        !> Index of q = 0
        integer, allocatable :: indices_zero_q(:)
        !> First q-vector for current rank for mpi distribution
        integer :: iq_start
        !> First q-vector for current rank for mpi distribution
        integer :: iq_end

        n_qvecs = size(n_gqvecs)

        ngq_max = maxval(n_gqvecs)

        ! Find indices of zero and finite q-vectors, respectively
        indices_finite_q = indices_finite_vectors(qvecs)
        indices_zero_q = indices_zero_vectors(qvecs)
        call terminate_if_false(size(indices_finite_q) > 0, &
                                    'setup_dielectric_matrix_from_file: Bad q-grid.&
                                    Bad q-grid. Only contains Gamma, or points approximately equal to zero')
        ! Read only on root process
        if (mpiglobal%is_root) then

            iq = indices_zero_q(1)

            write (strnum, '(I5.5)') iq
            fname = trim(eps0_dirname)//'/EPS0_Q'//trim(strnum)//'.OUT'
            call geteps0_zero_q(iq=iq, qvec=qvecs(:, iq), fname=fname,&
            &  iw=1, w=0._dp, eps0=this%body(:n_gqvecs(iq), :n_gqvecs(iq), iq),&
                    & eps0wg=this%wings, eps0hd=this%head)

            do i = 1, size(indices_finite_q)
                iq = indices_finite_q(i)

                write (strnum, '(I5.5)') iq

                fname = trim(eps0_dirname)//'/EPS0_Q'//trim(strnum)//'.OUT'

                iq = indices_finite_q(i)
                call geteps0_finite_q(iq=iq, qvec=qvecs(:, iq), fname=fname,&
                &  iw=1, w=0._dp, eps0=this%body(:n_gqvecs(iq), :n_gqvecs(iq), iq))
            end do

        end if

        ! Broadcast full dielectric matrix on all processes
        call xmpi_bcast(mpiglobal, this%head)
        call xmpi_bcast(mpiglobal, this%wings)
        call xmpi_bcast(mpiglobal, this%body)

    end subroutine

    !> Compares two dielectric matrices. Helpful for evaluating the
    !> quality of the additive screening approach, i.e. the DM of the total
    !> obtained via the additive approach and the explicit formalism can be compared.
    !> q-vectors and (G+q)-vectors are expected to be stored in QPOINTS_SCR.OUT
    !> and GQPOINTS/..., respectively.
    subroutine compare_dielectric_matrix(eps_1_dirname, eps_2_dirname, n_qvecs, tol)

        use putgeteps0, only: geteps0_finite_q, geteps0_zero_q
        use math_utils, only: all_close, diag, all_zero
        use constants, only: zzero
        use modmpi, only: mpiglobal, distribute_loop

        use exciting_mpi, only: xmpi_bcast
        use mod_kpointset, only: Gk_set, k_set, &
                                 delete_Gk_vectors, delete_k_vectors
        use expand_add_eps_helper, only: setup_gq_vecs_from_file, setup_q_vecs_from_file

        !> Directory where first matrix is stored
        character(*), intent(in) :: eps_1_dirname
        !> Directory where second matrix is stored
        character(*), intent(in) :: eps_2_dirname
        !> Number of q-vectors
        integer, intent(in) :: n_qvecs
        !> Tolerance
        real(dp), optional, intent(in) :: tol

        !> Number of (G+q)-vectors for each q-vector
        type(k_set) :: q_set
        !> (G+q)-vectors for each q-vector
        type(Gk_set) :: gq_set
        !> Dielectric matrix 1
        type(dielectric_matrix_type) :: eps_1
        !> Dielectric matrix 2
        type(dielectric_matrix_type) :: eps_2
        !> Filename
        character(256) :: fname1, fname2
        !> Full filename
        character(256) :: fname
        !> Directory  name
        character(256) :: dirname, syscommand
        !> Running index q-vectors
        integer :: iq, i

        !> Tolerance
        real(dp) :: tol_
        !> string for filename
        character(256) :: strnum
        !> Spin channel (needed to use 'Gk_set' type)
        integer, parameter :: ispin = 1

        tol_ = tol_default
        if (present(tol)) tol_ = tol

        ! Read q and (G+q)-vectors for heterostructure from file
        call setup_q_vecs_from_file('QPOINTS_SCR.OUT', n_qvecs, q_set)
        call setup_gq_vecs_from_file('GQPOINTS', 'QPOINTS_SCR.OUT', 1, n_qvecs, gq_set)
        call eps_1%init(gq_set%ngk(ispin, :), 1)
        call eps_2%init(gq_set%ngk(ispin, :), 1)

        call eps_1%setup_from_file(eps_1_dirname, &
                                    & gq_set%ngk(ispin, :), q_set%vkc)
        call eps_2%setup_from_file(eps_2_dirname,&
                                    & gq_set%ngk(ispin, :), q_set%vkc)
        if (mpiglobal%is_root) then
            write (*, *)
            write (*, *) '------Comparison dielectric matrices------'
            write (*, *) 'Directories:', trim(eps_1_dirname), '   ', trim(eps_2_dirname)
            write (*, *)

            ! Compare body for all q-vectors
            do iq = 1, n_qvecs

                write (*, *) all_close(eps_1%body(:gq_set%ngk(ispin, iq), :gq_set%ngk(ispin, iq), iq), &
                                       eps_2%body(:gq_set%ngk(ispin, iq), :gq_set%ngk(ispin, iq), iq), tol_)
            end do

            write (*, *)

            ! Compare head and wings for q=0
            write (*, *) all_close(eps_1%head, eps_2%head, tol_)
            write (*, *) all_close(eps_1%wings(:, 1, :), eps_2%wings(:, 1, :), tol_)
            write (*, *) all_close(eps_1%wings(:, 2, :), eps_2%wings(:, 2, :), tol_)

        end if

        call eps_1%delete
        call eps_2%delete
        call delete_Gk_vectors(gq_set)
        call delete_k_vectors(q_set)

    end subroutine


    !> A wrapper for puteps0_finite_q and  puteps0_zero_q:
    !> Writes the dielectric matrix for all q-vectors employed in a BSE
    !> calculation to binary files. Files read for zero frequency.
    !> The directory name has to be specified, while the files
    !> are named EPS0_QXXXXX.OUT, where XXXXX is the q-point index
    !> with leading zeros.
    subroutine write_dielectric_matrix(this, eps0_dirname, &
                                       n_gqvecs, qvecs)
        use modmpi, only: mpiglobal, terminate_if_false, distribute_loop
        use math_utils, only: all_zero
        use constants, only: zzero
        use putgeteps0, only: puteps0_finite_q, puteps0_zero_q
        use grid_utils, only: indices_zero_vectors, indices_finite_vectors
        !> Directory where dielectric matrices are stored
        character(*), intent(in) :: eps0_dirname
        !> Number of (G+q)-vectors for all q-vectors
        integer, intent(in) :: n_gqvecs(:)
        !> q-vectors
        real(dp), intent(in) :: qvecs(:, :)
        !> Dielectric matrix
        class(dielectric_matrix_type), intent(inout) :: this

        !> Number of q-vectors
        integer :: n_qvecs
        !> Running index (G+q)-vectors in unit cell
        integer :: iq, i
        !> Name of a file
        character(256) :: fname
        !> string for filename
        character(256) :: strnum
        !> Indices of finite q
        integer, allocatable :: indices_finite_q(:)
         !> Indices of finite q
        integer, allocatable :: indices_zero_q(:)
        !> Index of q = 0
        integer :: index_zero_q
        !> First q-vector for current rank for mpi distribution
        integer :: iq_start
        !> First q-vector for current rank for mpi distribution
        integer :: iq_end

        n_qvecs = size(n_gqvecs)

        ! Find indices of zero and finite q-vectors
        indices_finite_q = indices_finite_vectors(qvecs)
        indices_zero_q = indices_zero_vectors(qvecs)
        index_zero_q =indices_zero_q(1)
        call terminate_if_false(size(indices_finite_q) > 0, &
                                'write_dielectric_matrix: Bad q-grid.')

        ! Write for q = 0 only on root process
        if (mpiglobal%is_root) then

            iq = index_zero_q

            write (strnum, '(I5.5)') iq
            fname = trim(eps0_dirname)//'/EPS0_Q'//trim(strnum)//'.OUT'
            call puteps0_zero_q(iq=iq, qvec=qvecs(:, iq), fname=fname,&
            &  iw=1, w=0._dp, eps0=this%body(:n_gqvecs(iq), :n_gqvecs(iq), iq),&
                    & eps0wg=this%wings, eps0hd=this%head)

        end if

        call distribute_loop(mpiglobal, size(indices_finite_q), iq_start, iq_end)

        do i = iq_start, iq_end
            iq = indices_finite_q(i)

            write (strnum, '(I5.5)') iq

            fname = trim(eps0_dirname)//'/EPS0_Q'//trim(strnum)//'.OUT'

            iq = indices_finite_q(i)
            call puteps0_finite_q(iq=iq, qvec=qvecs(:, iq), fname=fname,&
            &  iw=1, w=0._dp, eps0=this%body(:n_gqvecs(iq), :n_gqvecs(iq), iq))
        end do
    end subroutine
end module