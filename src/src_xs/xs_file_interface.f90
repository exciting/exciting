!> Routines for reading objects written by an xs calculation.
!> TODO (Max / Benedikt): Check which routines similar to this one
!> will be merged with the fast-BSE code and remove duplicates.
module xs_file_interface
    use precision, only: sp, dp
    use m_getunit, only: getunit

    implicit none

    private

    public :: read_cartesian_coordinates, &
              read_n_gq, write_q_vectors,&
              read_gq_vectors, write_gq_vectors
contains

    !> Reads coordinates of the grid points from an exciting grid point
    !> file such as `KPOINTS.OUT` etc.
    subroutine read_cartesian_coordinates(file_name, cartesian_coordinates)
        !> File name
        character(*), intent(in) :: file_name

        real(dp), intent(out) :: cartesian_coordinates(:, :)

        integer :: i, n_vecs
        integer :: dummy_int
        integer :: unit_id
        real(dp) :: dummy_vec(3)

        character(256) :: charline

        call getunit(unit_id)

        open (unit=unit_id, file=trim(file_name), status='old', action='read')
        read (unit_id, *) n_vecs

        do i = 1, n_vecs
            read (unit_id, *) dummy_int, dummy_vec, cartesian_coordinates(:, i)
        end do

        close (unit_id)
    end subroutine read_cartesian_coordinates

    !> Extract the number of \(|\mathbf{G}+\mathbf{q}|\) points
    !> per \(\mathbf{q}\)-point from file and distributes it to all processes.
    subroutine read_n_gq(file_name_q, n_gqvecs)

        use modmpi, only: mpiglobal
        use exciting_mpi, only: xmpi_bcast

        !> Name of file containing q-vectors
        character(*), intent(in) :: file_name_q
        !> Number of q-vectors
        integer(sp), intent(out) :: n_gqvecs(:)

        !> Running index q-vectors
        integer :: i_q
        !> Number of q-vectors
        integer ::  n_qvecs
        !> Integer for reading
        integer :: dummy_int
        !>  Skipped vector for reading
        real(dp) :: vec(6)
        integer :: unit_id

        call getunit(unit_id)

        open (unit=unit_id, file=trim(file_name_q), status='old', action='read')
        read (unit_id, *) n_qvecs
        do i_q = 1, n_qvecs
            read (unit_id, *) dummy_int, vec, n_gqvecs(i_q)

        end do
        close (unit_id)

    end subroutine read_n_gq

    !> Wrapper for read_n_gq and read_latt_coords: Reads
    !> all q-vectors and all (G+q)-vectors from file and broadcasts them
    !> to all MPI ranks. The name of the file containing the q-vectors and
    !> the name of the directory with the files containing (G+q)-vectors
    !> has to be given. The names of the files containing the (G+q)-vectors
    !> is assumed to be GQPOINTS_SCR_QXXXXX.OUT where XXXXX is the q-vector index
    !> with leading zeros.
    subroutine read_gq_vectors(dir_name_Gq, n_gqvecs, gq_vecs)
        use modmpi, only: mpiglobal, mpi_allgatherv_ifc
        use exciting_mpi, only: xmpi_bcast
        use modxs, only: qpari, qparf

        !> Directory with the files containing the (G+q)-vectors
        character(*), intent(in) :: dir_name_Gq
        !> Number of (G+q)-vectors  for each q-vector
        integer, intent(in) :: n_gqvecs(:)

        !>  (G+q)-vectors for each q-vector
        real(dp), intent(out) :: gq_vecs(:, :, :)
        !> Number of q-vectors
        integer :: n_qvecs
        !> string for filename
        character(256) :: strnum
        !> Filename of file containing (G+q)-vectors
        character(256) :: fname
        !> Running index q-vectors
        integer :: iq
        !> Maximum number of (G+q)-vectors
        integer :: ngq_max

        n_qvecs = size(n_gqvecs)
        ngq_max = maxval(n_gqvecs)

        n_qvecs = size(n_gqvecs)

        do iq = 1, n_qvecs
            write (strnum, '(I5.5)') iq
            fname = trim(dir_name_Gq)//'/'//'GQPOINTS_SCR_Q'//trim(strnum)//'.OUT'
            call read_cartesian_coordinates(fname, gq_vecs(:, :n_gqvecs(iq), iq))

        end do

    end subroutine


    !> Writes the (G+q)-vectors for one given q-vector to file.
    !> The directory name has to be specified as input, while the file name
    !> is set to be GQPOINTS_SCR_QXXXXX.OUT, where XXXXX is the q-vector index
    !> with leading zeros.
    subroutine write_gq_vectors(gqvecs_dirname, iq, gq_vecs_lat, &
                                gq_vecs_cart, n_gq)

        use m_getunit, only: getunit
        !> Directory name where (G+q)-vectors are stored
        character(*), intent(in) :: gqvecs_dirname
        !> Running index q-vector
        integer :: iq
        !> (G+q)-vectors  in lattice coordinates
        real(dp), intent(in) :: gq_vecs_lat(:, :)
        !> (G+q)-vectors in cartesian coordinates
        real(dp), intent(in) :: gq_vecs_cart(:, :)
        !> File id
        integer :: fid
        !> Status for I/O
        integer :: stat
        !> Running index (G+q)-vector
        integer(sp) :: igq
        !> File name
        character(256) :: fname
        !> Strin for q-vector index
        character(256) :: strnum
        !> Number of (G+q)-vectors
        integer(sp) :: n_gq

        ! Generate filename and open file
        write (strnum, '(I5.5)') iq
        fname = trim(gqvecs_dirname)//'/GQPOINTS_SCR_Q'//trim(strnum)//'.OUT'
        call getunit(fid)
        open (unit=fid, file=trim(fname), action='write', iostat=stat)

        ! Write header
        write (fid, '(I6, " : ngq; G+q-point, vgql, vgqc, gqc, |G| below")') n_gq

        ! Write (G+q)-vectors
        do igq = 1, n_gq
            write (fid, '(I6, 3G18.10,2x,3G18.10,2x,2G18.10)') &
                igq, gq_vecs_lat(:, igq), gq_vecs_cart(:, igq), &
                norm2(gq_vecs_cart(:, igq)), &
                norm2(gq_vecs_cart(:, igq) - gq_vecs_cart(:, 1))
        end do

        close (fid)
    end subroutine

    !> Writes all q-vectors to file. The filename has to be specified as input.
    subroutine write_q_vectors(qvecs_fname, qvecs_lat, qvecs_cart, n_gqvecs)
        use m_getunit, only: getunit
        use modmpi, only: terminate_mpi_env, mpiglobal
        integer :: fid
        !> q-vectors in lattice coordinates
        real(dp), intent(in) :: qvecs_lat(:, :)
        !> q-vectors in cartesian coordinates
        real(dp), intent(in) :: qvecs_cart(:, :)
        !> Number of (G+q)-vectors per q-vector
        integer, intent(in) :: n_gqvecs(:)
        !> Filename for q-vectors
        character(*), intent(in) :: qvecs_fname
        !> I/O status
        integer(sp) :: stat
        !> Running index q-vector
        integer :: iq
        !> Number of q-vectors
        integer :: n_qvecs

        n_qvecs = size(n_gqvecs)

        ! Open file
        call getunit(fid)
        open (unit=fid, file=trim(qvecs_fname), action='write', iostat=stat)
        if (stat /= 0) then
            write (*, *) "Error opening file, iostat:", stat, qvecs_fname
            call terminate_mpi_env(mpiglobal)
        end if

        ! Write header
        write (fid, '(I6, " : nqpt; q-point, vql, vqc, ngq below")') n_qvecs

        ! Write q-vectors
        do iq = 1, n_qvecs
            write (fid, '(i6, 6E18.10, i8)') iq, qvecs_lat(:, iq), &
                qvecs_cart(:, iq), n_gqvecs(iq)
        end do

        close (fid)
    end subroutine write_q_vectors
end module
