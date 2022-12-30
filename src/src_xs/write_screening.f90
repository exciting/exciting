!> Module for writing the dielectric matrix (DM),
!>
!> \[
!>  \varepsilon_{\mathbf{G}\mathbf{G}'}(\mathbf{q}),
!>                                                  \]
!> and the screened Coulomb interaction (SCI),
!>
!> \[
!>  W_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) =
!>    \varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q})
!>    v_{\mathbf{G}'}(\mathbf{q}),
!>                                                       \]
!> on the full set of non-reduced q-vectors.
!> In both cases, the starting point is the DM on the set of
!> reduced q-vectors, which can be related to the DM on the non-reduced
!> q-grid by the use of symmetry operations.
!> The DM is written to direct access binary files,
!> while the SCI is written to HDF5.
!>
!> Adapted from the procedure used in scrcoulint.f90.
module write_screening

    use constants, only: zzero
    use precision, only: i32, dp

    implicit none

    private

    public :: write_screening_launcher

contains

    !> Launcher routine for the generation of the dielectric matrix (DM)
    !> and the screened Coulomb interaction on the
    !> full set of non-reduced q-vectors. This routine should be called
    !> after computing the DM on the set of reduced q-vectors (task="screen"
    !> with xs-flag reduceq="true").
    !>
    !> The following steps are done in both cases:
    !>  i) Initialization of reduced and non-reduced q-grids and the
    !>      corresponding (G+q)-grids. The q-grids and (G+q)-grids are
    !>      stored in the global modxs variables 'q' and 'g_q', respectively.
    !>
    !>  ii) Writing the non-reduced q-grid and the corresponding
    !>       (G+q)-grid to file (to file QPOINTS_SCR.OUT
    !>       and directory GQPOINTS/, respectively).
    !>
    !> The following step is done only if the task is to write the DM:
    !>
    !>  iii) Reading of the DM on the reduced q-grid and passing it to the
    !>       routine 'write_dielectric_matrix'.
    !>
    !> To be consistent with the naming conventions chosen for the types
    !> 'gk_set' and 'q_set', the suffix '_nr' is assigned to all variables
    !> referring to quantities on the non-reduced q-grid.
    subroutine write_screening_launcher(task)
        use mod_xsgrids, only: xsgrids_init, q, g_q
        use mod_Gkvector, only: gkmax
        use modmpi, only: terminate_mpi_env, mpiglobal
        use modinput, only: input
        use modxs, only: qpari, qparf, ivqr, iqmapr, &
                         vqlr, vqcr, wqptr, ngqr, &
                         eps0dirname, nqptr
        use mod_qpoint, only: iqmap, vql, vqc, ivq, nqpt
        use xs_file_interface, only: write_gq_vectors, write_q_vectors

        !> Defines if DM or screened Coulomb interaction should be written
        character(*), intent(in) :: task
        !> Directory name for (G+q)-vectors for non-reduced q-vectors
        character(256), parameter :: gqvecs_dirname = 'GQPOINTS'
        !> File name for non-reduced q-vectors
        character(256), parameter :: qvecs_fname = 'QPOINTS_SCR.OUT'
        !>  Folder name for dielectric matrix for non-reduced q-vectors
        character(256), parameter :: eps0_dirname = 'EPS0'

        !> Gamma point
        real(dp) :: gamma_point(3) = [0._dp, 0._dp, 0._dp]
        !> Spin index, only needed to use types q_set and gk_set
        integer(i32), parameter :: ispin = 1
        !> Running index non-reduced q-vectors
        integer(i32) :: iq_nr
        !> Number of (G+q)-vectors for one non-reduced q-vector
        integer(i32) :: ngq_nr
        !> System command
        character(256) :: syscommand
        !> Error message
        character(256) :: error_message

        !> Head of DM
        complex(dp), allocatable  :: eps_head(:, :)
        !> Wings of DM
        complex(dp), allocatable  :: eps_wings(:, :, :)
        !> Body of DM for all reduced q-vectors
        complex(dp), allocatable  :: eps_body(:, :, :)
        !> True if the inquired file exists
        logical :: file_exists

        ! Initialize q-grid and (G+q)-grid for both the reduced and the
        ! non-reduced set of q-vectors
        gamma_point = (/0._dp, 0._dp, 0._dp/)

        call init0
        call init1
        call xssave0
        call init2
        call xsgrids_init(gamma_point, gkmax, makegq_=.true.)
        call init1offs(gamma_point)
        call init2offs(gamma_point, input%xs%reduceq)

        ! Copy q-point results from mod_qpoint into modxs variables
        nqptr = nqpt
        ivqr = ivq
        iqmapr = iqmap
        vqlr = vql
        vqcr = vqc
        call init2offs(gamma_point, .false.)

        ! Write q-grid and (G+q)-grid to file for all non-reduced q-vectors
        if (mpiglobal%is_root) then

            call write_q_vectors(qvecs_fname, q%qset%vklnr, &
                                 q%qset%vkcnr, g_q%ngknr(ispin, :))

            do iq_nr = 1, q%qset%nkptnr

                ngq_nr = g_q%ngknr(ispin, iq_nr)

                call write_gq_vectors(gqvecs_dirname, iq_nr, &
                                      g_q%vgknrl(:, :, ispin, iq_nr), &
                                      g_q%vgknrc(:, :, ispin, iq_nr), ngq_nr)
            end do
        end if

        select case (trim(task))

        case ('write_dielectric_matrix')

            ! Read dielectric matrix on reduced q-vectors
            call get_dielectric_matrix(trim(adjustl(eps0_dirname)), &
                                       g_q%ngk(ispin, :), q%qset%vkc, &
                                       eps_head, eps_wings, eps_body)

            call write_dielectric_matrix(eps_head, eps_wings, &
                                         eps_body)

        case ('write_screened_coulomb')

            call write_screened_coulomb_interaction

        case default
            if (mpiglobal%is_root) then
                error_message = 'write_screening_init: Invalid task: '//trim(task)
                call terminate_mpi_env(mpiglobal, error_message)
            end if
        end select
    end subroutine

    !> Given the dielectric matrix on the set of symmetry-reduced q-vectors,
    !> this routine generates the matrix on all q-vectors and writes it to files
    !> into the directory EPS0/EPS0_QXXXXX.OUT,
    !> where XXXXX is the q-point indec with leading zeros.
    !> The matrix on the full q-grid is obtained by
    !> making use of the symmetry relations described in the module
    !> genphasedm.f90.
    !> To be consistent with the naming conventions chosen for the types
    !> 'gk_set' and 'q_set', the suffix '_nr' is assigned to all variables
    !> referring to quantities on the non-reduced q-grid.
    subroutine write_dielectric_matrix(eps_head, eps_wings, &
                                       eps_body)

        use mod_xsgrids, only: q, g_q
        use modxs, only: qpari, qparf
        use mod_symmetry, only: maxsymcrys
        use modinput, only: input
        use m_getunit, only: getunit
        use putgeteps0, only: puteps0_finite_q, puteps0_zero_q
        use modmpi, only: mpiglobal, terminate_mpi_env

        !> Head of DM at Gamma for reduced q-vectors
        complex(dp), intent(in)  :: eps_head(:, :)
        !> Wings of DM for reduced q-vectors
        complex(dp), intent(in)  :: eps_wings(:, :, :)
        !> Body of reduced DM for reduced q-vectors
        complex(dp), intent(in) :: eps_body(:, :, :)

        !> Running index reduced q-vectors
        integer(i32) :: iq
        !> Running index non-reduced q-vectors
        integer(i32) :: iq_nr
        !> Number of symmetries transforming non-reduced q-point to reduced one
        integer(i32) :: nsc, ivgsym(3)
        !> Symmtetry maps
        integer(i32) :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
        !> Index symmetry operation
        integer(i32) :: jsym
        !> Index symmetry operation
        integer(i32) :: jsymi
        !> Map between (G+q)-vectors of reduced and non-reduced q-grid
        integer(i32), allocatable :: igqmap(:)
        !> Phasefactors for dielectric matrix
        complex(8), allocatable :: phasefactors(:, :)
        !> Checks if non-trivial phase appears at least for one (G,Gp) component
        logical :: tphf
        !> Body of DM for one non-reduced q-vector
        complex(dp), allocatable  :: eps_body_nr(:, :)
        !> Filename
        character(256) :: fname
        !> string for filename
        character(256) :: strnum
        !> Indices of finite q
        integer(i32), allocatable :: indices_finite_q(:)
        !> Index of q = 0
        integer(i32), allocatable :: indices_zero_q(:)
        !> Dummy index for creating integer array
        integer(i32) :: i
        !> One non-reduced q-vector
        real(dp) :: q_vec_nr(3)
        !> One reduced q-vector
        real(dp) :: q_vec(3)
        !> Number of (G+q)-vectors for one non-reduced q-vector
        integer(i32) :: ngq_nr
        !> Spin index, only needed to use types q_set and gk_set
        integer(i32), parameter :: ispin = 1
        !>  Directory name for dielectric matrix for non-reduced q-vectors
        character(256), parameter :: eps0_dirname = 'EPS0'

        ! Find indices of zero and finite q-vectors, respectively
        indices_finite_q = pack([(i, i=1, q%qset%nkptnr)], norm2(q%qset%vkcnr, dim=1) > 1.e-9_dp)
        indices_zero_q = pack([(i, i=1, q%qset%nkptnr)], norm2(q%qset%vkcnr, dim=1) < 1.e-9_dp)
        if (size(indices_zero_q) > 1) then
            call terminate_mpi_env(mpiglobal, &
                                   'Error eps_on_all_q: More than one &
                                   non-reduced q-vector is zero.')
        end if

        ! Only need to write q=0 from one process. No symmetry operation needed.
        if (mpiglobal%is_root) then

            ! Index and lattice coordinates of current non-reduced q-vector
            iq_nr = indices_zero_q(1)
            q_vec_nr = q%qset%vklnr(:, iq_nr)

            ! Number of (G+q)-vectors for non-reduced q-vector
            ngq_nr = g_q%ngknr(ispin, iq_nr)

            write (strnum, '(I5.5)') iq_nr
            fname = trim(adjustl(eps0_dirname))//'/EPS0_Q'//trim(strnum)//'.OUT'

            call puteps0_zero_q(iq=iq_nr, qvec=q_vec_nr, fname=fname,&
                    &  iw=1, w=0._dp, eps0=eps_body(:ngq_nr, :ngq_nr, iq_nr),&
                    & eps0wg=eps_wings, eps0hd=eps_head)
        end if

        ! Find results for finite non-reduced q-vectors with the help
        ! of the results for the corresponding reduced q-vector using
        ! symmetry operations.
        allocate (phasefactors(g_q%ngknrmax, g_q%ngknrmax))
        allocate (igqmap(g_q%ngknrmax), source=0)

        ! Distribute finite q-vectors of reduced q-grid over all processes
        ! Respective first / last index stored in qpari / qparf
        call genparidxran('q', size(indices_finite_q))

        do i = qpari, qparf

            ! Index and lattice coordinates of current non-reduced q-vector
            iq_nr = indices_finite_q(i)
            q_vec_nr = q%qset%vklnr(:, iq_nr)

            ! Index and lattice coordinates of reduced q-vector
            iq = q%qset%ik2ikp(iq_nr)
            q_vec = q%qset%vkl(:, iq)

            ! Number of (G+q)-vector for non-reduced q-vector
            ngq_nr = g_q%ngknr(ispin, iq_nr)

            ! Find a crystal symmetry operation that rotates the non-reduced
            ! (G+q_nr)-vectors  onto (G'+q)-vectors (where q is a
            ! vector from the reduced q-grid) and generate a Map G' --> G
            call findsymeqiv(input%xs%bse%fbzq, q_vec_nr, q_vec, nsc, sc, ivgsc)
            call findgqmap(iq_nr, iq, nsc, sc, ivgsc, ngq_nr, jsym, jsymi, ivgsym, igqmap(:ngq_nr))

            ! Generate phase factor for dielectric matrix due to non-primitive
            ! translations
            call genphasedm(iq_nr, jsym, g_q%ngknrmax, ngq_nr, &
                            phasefactors, tphf)

            ! Obtain dielectric matrix for non-reduced q-vector
            if (allocated(eps_body_nr)) deallocate (eps_body_nr)
            allocate (eps_body_nr(ngq_nr, ngq_nr))

            eps_body_nr = &
                phasefactors(:ngq_nr, :ngq_nr)* &
                eps_body(igqmap(:ngq_nr), igqmap(:ngq_nr), iq)

            ! Write to file
            write (strnum, '(I5.5)') iq_nr
            fname = trim(adjustl(eps0_dirname))//'/EPS0_Q'//trim(strnum)//'.OUT'

            call puteps0_finite_q(iq=iq_nr, qvec=q%qset%vkcnr(:, iq_nr), &
                                  fname=fname, iw=1, w=0._dp, eps0=eps_body_nr)

        end do

    end subroutine

    !> Generates the matrix of the screened Coulomb interaction
    !> \[
    !>  W_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) =
    !>    \varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q})
    !>    v_{\mathbf{G}'}(\mathbf{q}),
    !>                                                       \]
    !> for the full set of non-reduced q-vectors and writes it to
    !> the HDF5 file 'bse_output.h5'.
    !> The starting point is the DM on the set of reduced q-vectors,
    !> which can be related to the DM on the non-reduced
    !> q-grid by the use of symmetry relations as explained in
    !> the module genphasedm.f90.
    !> To be consistent with the naming convecntions chosen for the types
    !> 'gk_set' and 'q_set', the suffix '_nr' is assigned to all variables
    !> referring to quantities on the non-reduced q-grid.
    ! TODO (Max / Benedikt): It should be checked if this is still needed
    ! or if Benedikt has already implemented a newer version
    subroutine write_screened_coulomb_interaction

        use mod_xsgrids, only: q, g_q
        use mod_symmetry, only: maxsymcrys
        use modinput, only: input
        use m_getunit, only: getunit
        use modmpi, only: mpiglobal, terminate_mpi_env, mpi_allgatherv_ifc
        use modxs, only: eps0dirname, qpari, qparf
        use mod_hdf5, only: fhdf5, hdf5_exist_group, hdf5_create_group, &
                            hdf5_write, hdf5_delete_object

        !> Running index reduced q-vectors
        integer(i32) :: iq
        !> Running index non-reduced q-vectors
        integer(i32) :: iq_nr
        !> Number of symmetries transforming non-reduced q-point to reduced one
        integer(i32) :: nsc, ivgsym(3)
        !> Symmtetry maps
        integer(i32) :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
        !> Index symmtry operation
        integer(i32) :: jsym
        !> Index symmetry operation
        integer(i32) :: jsymi
        !> Map between (G+q)-vectors of reduced and non-reduced q-grid
        integer(i32), allocatable :: igqmap(:)
        !> Phasefactors for dielectric matrix
        complex(8), allocatable :: phasefactors(:, :)
        !> Checks if non-trivial phase appears at least for one (G,Gp) component
        logical :: tphf
        !> Screened Coulomb interaction (W_GG'(q)) for all reduced q-vectors
        complex(dp), allocatable  :: screened_coulomb(:, :, :)
        !> Screened Coulomb interaction (W_GG'(q)) for one non-reduced q-vector
        complex(dp), allocatable  :: screened_coulomb_nr(:, :, :)
        !> Filename
        character(256) :: fname
        !> string for filename
        character(256) :: strnum
        !> Indices of finite q
        integer(i32), allocatable :: indices_finite_q(:)
        !> Index of q = 0
        integer(i32), allocatable :: indices_zero_q(:)
        !> Dummy index for creating integer array
        integer(i32) :: i
        !> One non-reduced q-vector
        real(dp) :: q_vec_nr(3)
        !> One reduced q-vector
        real(dp) :: q_vec(3)
        !> Number of (G+q)-vectors for one non-reduced q-vector
        integer(i32) :: n_gq_nr
        !> Number of (G+q)-vectors for one reduced q-vector
        integer(i32) :: n_gq
        !> Spin index, needed to use types q_set and gk_set
        integer(i32), parameter :: ispin = 1
        ! HDF5 variables
        character(256) :: ciq, gname, group

        ! Set name for global variable, needed in gensclieff.f90
        eps0dirname = 'EPS0'

        allocate (igqmap(g_q%ngknrmax), source=0)
        allocate (screened_coulomb(g_q%ngkmax, g_q%ngkmax, q%qset%nkpt))

        ! Distribute reduced q-grid over all processes
        ! Respective first / last index stored in qpari / qparf
        call genparidxran('q', q%qset%nkpt)

        ! Compute screened Coulomb interaction for reduced q-grid
        do iq = qpari, qparf
            n_gq = g_q%ngk(ispin, iq)
            call genscclieff(iq, g_q%ngkmax, n_gq, screened_coulomb(:, :, iq))
        end do

        ! Communicate array-parts wrt. reduced q-grid
        call mpi_allgatherv_ifc(set=q%qset%nkpt, rlen=g_q%ngkmax**2,&
          & zbuf=screened_coulomb, inplace=.true., comm=mpiglobal)

        ! Find results for finite non-reduced q-vectors
        allocate (phasefactors(g_q%ngkmax, g_q%ngkmax))

        allocate (screened_coulomb_nr(g_q%ngkmax, g_q%ngkmax, q%qset%nkptnr))

        ! Distribute non-reduced q-grid over all processes
        ! Respective first / last index stored in qpari / qparf
        call genparidxran('q', q%qset%nkptnr)

        ! Obtain screened Coulomb interaction for non-reduced q-grid
        ! with the help of the results for the corresponding
        ! reduced q-vector using  symmetry operations.
        do iq_nr = qpari, qparf

            ! Lattice coordinates of current non-reduced q-vector
            q_vec_nr = q%qset%vklnr(:, iq_nr)

            ! Index and lattice coordinates of reduced q-vector
            iq = q%qset%ik2ikp(iq_nr)
            q_vec = q%qset%vkl(:, iq)

            ! Number of (G+q)-vector for non-reduced q-vector
            n_gq_nr = g_q%ngknr(ispin, iq_nr)

            ! Find a crystal symmetry operation that rotates the non-reduced
            ! (G+q_nr)-vectors  onto (G'+q)-vectors (where q is a
            ! vector from the reduced q-grid) and generate a Map G' --> G
            call findsymeqiv(input%xs%bse%fbzq, q_vec_nr, q_vec, nsc, sc, ivgsc)
            call findgqmap(iq_nr, iq, nsc, sc, ivgsc, n_gq_nr, jsym, jsymi, ivgsym, igqmap(:n_gq_nr))

            ! Generate phase factor for dielectric matrix due to non-primitive
            ! translations
            call genphasedm(iq_nr, jsym, g_q%ngknrmax, n_gq_nr, &
                            phasefactors, tphf)

            screened_coulomb_nr(:n_gq_nr, :n_gq_nr, iq_nr) = &
                phasefactors(:n_gq_nr, :n_gq_nr)* &
                screened_coulomb(igqmap(:n_gq_nr), igqmap(:n_gq_nr), iq)

        end do

        ! Communicate array-parts wrt. reduced q-grid
        call mpi_allgatherv_ifc(set=q%qset%nkptnr, rlen=g_q%ngkmax**2,&
          & zbuf=screened_coulomb_nr, inplace=.true., comm=mpiglobal)

        ! Write to HDF5
#ifdef _HDF5_
        if (mpiglobal%is_root) then
            if (.not. hdf5_exist_group(fhdf5, '/', 'screenedpotential')) then
                call hdf5_create_group(fhdf5, '/', 'screenedpotential')
            end if
            gname = "/screenedpotential"

            if (hdf5_exist_group(fhdf5, '/', gname)) call hdf5_delete_object(fhdf5, '/', gname)
            call hdf5_create_group(fhdf5, '/', gname)

            ! loop over all non-reduced q-vectors
            do iq_nr = 1, q%qset%nkptnr
                write (ciq, '(I4.4)') iq_nr
                if (.not. hdf5_exist_group(fhdf5, trim(adjustl(gname)), ciq)) then
                    call hdf5_create_group(fhdf5, trim(adjustl(gname)), ciq)
                end if
                group = "/screenedpotential/"//ciq//'/'
                call hdf5_write(fhdf5, group, "wqq", screened_coulomb_nr(1, 1, iq_nr), shape(screened_coulomb_nr(:, :, iq_nr)))
            end do
        end if
#endif

    end subroutine

    !> Wrapper for geteps0_finite_q and geteps0_zero_q:
    !> Reads dielectric matrix for all q-vectors from file.
    subroutine get_dielectric_matrix(eps0_dirname, n_gqvecs, qvecs_cart, &
                                     eps_head, eps_wings, eps_body)
        use math_utils, only: all_zero
        use constants, only: zzero
        use putgeteps0, only: geteps0_finite_q, geteps0_zero_q
        use modmpi, only: mpiglobal, terminate_mpi_env, mpi_allgatherv_ifc
        use modxs, only: qpari, qparf
        use exciting_mpi, only: xmpi_bcast

        !> Directory where dielectric matrices are stored
        character(*), intent(in) :: eps0_dirname
        !> Number of (G+q)-vectors per q-vectors
        integer(i32), intent(in) :: n_gqvecs(:)
        !> q-vectors in cartesian coordinates
        real(dp), intent(in) :: qvecs_cart(:, :)
        !> Head of dielectric matrix at Gamma
        complex(dp), allocatable, intent(out) :: eps_head(:, :)
        !> Wings of dielectric matrix at Gamma
        complex(dp), allocatable, intent(out) :: eps_wings(:, :, :)
        !> Body of dielectric matrix at all q-vectors
        complex(dp), allocatable, intent(out)  :: eps_body(:, :, :)

        !> Number of q-vectors
        integer(i32) :: n_qvecs
        !> Running index (G+q)-vectors in unit cell
        integer(i32) :: iq
        !> Name of a file
        character(256) :: fname
        !> string for filename
        character(256) :: strnum
        !> Maximum number of (G+q)-vectors per q-vector
        integer(i32) :: n_gq_max
        !> Indices of finite q
        integer(i32), allocatable :: indices_finite_q(:)
        !> Index of q = 0
        integer(i32), allocatable :: indices_zero_q(:)
        !> Dummy index for creating integer array
        integer(i32) :: i

        n_gq_max = maxval(n_gqvecs)

        n_qvecs = size(n_gqvecs)

        ! Find indices of zero and finite q-vectors, respectively
        indices_finite_q = pack([(i, i=1, n_qvecs)], norm2(qvecs_cart, dim=1) > 1.e-9_dp)
        indices_zero_q = pack([(i, i=1, n_qvecs)], norm2(qvecs_cart, dim=1) < 1.e-9_dp)

        allocate (eps_body(n_gq_max, n_gq_max, n_qvecs), source=zzero)
        allocate (eps_head(3, 3))
        allocate (eps_wings(n_gqvecs(indices_zero_q(1)), 2, 3))

        if (mpiglobal%is_root) then

            ! Read DM for q = 0
            iq = indices_zero_q(1)
            write (strnum, '(I5.5)') iq
            fname = trim(eps0_dirname)//'/EPS0_Q'//trim(strnum)//'.OUT'

            call geteps0_zero_q(iq=iq, qvec=qvecs_cart(:, iq), fname=fname,&
                &  iw=1, w=0._dp, eps0=eps_body(:n_gqvecs(iq), :n_gqvecs(iq), iq),&
                        & eps0wg=eps_wings, eps0hd=eps_head)

            ! Read DM for q /= 0
            do i = 1, size(indices_finite_q)
                iq = indices_finite_q(i)
                write (strnum, '(I5.5)') iq
                fname = trim(eps0_dirname)//'/EPS0_Q'//trim(strnum)//'.OUT'

                call geteps0_finite_q(iq=iq, qvec=qvecs_cart(:, iq), fname=fname,&
                 &  iw=1, w=0._dp, eps0=eps_body(:n_gqvecs(iq), :n_gqvecs(iq), iq))
            end do
        end if

        ! Broadcast dielectric matrix to all processes
        call xmpi_bcast(mpiglobal, eps_head)
        call xmpi_bcast(mpiglobal, eps_wings)
        call xmpi_bcast(mpiglobal, eps_body)

    end subroutine

    

end module write_screening
