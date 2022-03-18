module phonon_screening_file_interface

    use precision, only: sp, dp
    use modmpi, only: mpiglobal, terminate_mpi_env

    implicit none

contains

    !> Read dielectric tensor and Born effective charges
    !> from Quantum Espresso ph.x output (xxx.dyn1)
    !> xxx.dyn1 is the dynamical matrix at the Gamma point
    subroutine qe_read_eps_inf_zstar(natoms, fname, eps_infty, zstar)
        use m_getunit, only: getunit
        use errors_warnings, only: terminate_if_false
        use modmpi, only: mpiinfo

        ! I/O
        !> number of atoms
        integer(sp), intent(in) :: natoms(:)
        !> filename
        character(len=40), intent(in) :: fname
        !> High frequency dielectric tensor \varepsilon_\infty
        real(dp), intent(out)   ::  eps_infty(:, :)
        !> Born effective charge tensor
        real(dp), intent(out)   :: zstar(:, :, :)

        ! Help variables
        type(mpiinfo) :: mpiglobal
        !> Number of species
        integer(sp) :: nspecies
        !> Total number of atoms
        integer(sp) :: natmtot
        !> Running index lines
        integer(sp) :: ilines
        !> fileindex
        integer(sp) :: fid
        !> Running index atoms
        integer(sp)  :: iatom
        !> Running index species
        integer(sp)  :: ispecies
        !> Running index combined atoms and species
        integer(sp)  :: ias

        ! Get total number of atoms
        natmtot = size(zstar, dim=3)

        ! Get number of species
        nspecies = size(natoms)

        call getunit(fid)
        open (unit=fid, file=trim(fname), status='old', action='read')

        ! Skip file infos
        do ilines = 1, 7
            read (fid, *)
        end do

        ! Skip species info
        do ispecies = 1, nspecies
            read (fid, *)
        end do

        ! Skip atom info
        do ispecies = 1, nspecies
            do iatom = 1, natoms(ispecies)
                read (fid, *)
            end do
        end do

        ! Skip dynamical matrix info
        do ilines = 1, 5
            read (fid, *)
        end do

        ! Skip dynamical matrix
        do ilines = 1, 4*natmtot**2
            read (fid, *)
        end do

        ! Skip whitespace
        do ilines = 1, 3
            read (fid, *)
        end do

        ! Read dielectric tensor
        do ilines = 1, 3
            read (fid, *) eps_infty(ilines, :)
        end do

        ! Skip whitespace
        do ilines = 1, 3
            read (fid, *)
        end do

        ias = 0
        ! Read Born effective charges
        do ispecies = 1, nspecies
            do iatom = 1, natoms(ispecies)

                ias = ias + 1
                !Skip atom info
                read (fid, *)

                do ilines = 1, 3
                    read (fid, *) zstar(ilines, :, ias)
                end do
            end do
        end do

        close (fid)

        ! Check if number of atoms was correct
        call terminate_if_false(mpiglobal, ias == natmtot, &
                                error_message=' Total number of atoms not &
                                correct when reading: '//trim(fname))

    end subroutine qe_read_eps_inf_zstar

    !> Read phonon frequencies and eigenvectors
    !> from Quantum Espressomatdyn.x output (xxx.eig)
    !> xxxx.eig is the dynamical matrix at the Gamma point
    subroutine qe_read_phonon(natoms, fname, freq_ph, evec_ph_out)

        use m_getunit, only: getunit
        use unit_conversion, only: hartree_to_thz
        use errors_warnings, only: terminate_if_false
        use modmpi, only: mpiinfo

        ! I/O
        !> Number of atoms per species
        integer(sp), intent(in) :: natoms(:)
        !> filename
        character(len=40), intent(in) :: fname
        !> phonon frequencies
        real(dp), intent(out)  :: freq_ph(:, :)
        !> phonon eigenvectors for output
        complex(dp), intent(out), optional  :: evec_ph_out(:, :, :, :)
        !> phonon eigenvectors for output
        complex(dp), allocatable  :: evec_ph_local(:, :, :, :)
        !Help variables
        !> Running index qpoint
        integer(sp) :: iq
        !> Number of qpoints
        integer(sp) :: nqpoints
        !> Total number of atoms
        integer(sp) :: natmtot
        !> Number of species
        integer(sp) :: nspecies
        !> Running index phonon mode
        integer(sp) ::  imode
        !> Running index lines
        integer(sp) :: ilines
        !> Running index atoms
        integer(sp)  :: iatom
        !> Running index species
        integer(sp)  :: ispecies
        !> Running index combined atoms and species
        integer(sp)  :: ias
        !> fileindex
        integer(sp) :: fid
        !> Characters to read
        character(len=40) :: char(4)
        !> helper array phonon eigenvectors
        real(dp) :: e_real(3)
        !> helper array phonon eigenvectors
        real(dp) :: e_imag(3)
        integer(sp) :: var1
        !> Format for reading
        character(len=40) ::fmt
        !> mpiinfo
        type(mpiinfo) :: mpiglobal
        !> Number of phonon modes
        integer(sp) :: n_phonon_modes

        ! Get some infos
        natmtot = size(freq_ph, dim=1)/3._dp
        nspecies = size(natoms)
        nqpoints = size(freq_ph, dim=2)
        n_phonon_modes = 3*natmtot

        allocate (evec_ph_local(3, natmtot, n_phonon_modes, nqpoints))
        ! Define format for reading QE output file

        if (natmtot .le. 3._dp) then
            fmt = '(4X,A4,X,A1,4X,I1,A1,A1,7X,F10.5)'
        else
            fmt = '(4X,A4,X,A1,4X,I2,A1,A1,7X,F10.5)'
        end if

        call getunit(fid)
        open (unit=fid, file=trim(fname), status='old', action='read')

        do iq = 1, nqpoints
            ! Skip file info
            do ilines = 1, 4
                read (fid, *)
            end do

            !Read frequencies
            do imode = 1, n_phonon_modes

                read (fid, fmt) char(1), char(2), var1, char(3), char(4), freq_ph(imode, iq)
                ias = 0
                do ispecies = 1, nspecies
                    do iatom = 1, natoms(ispecies)
                        ias = ias + 1

                        ! read eigenvectors
                        read (fid, *) char(1), e_real(1), e_imag(1), e_real(2), e_imag(2), e_real(3), e_imag(3)

                        evec_ph_local(:, ias, imode, iq) = cmplx(e_real, e_imag)
                    end do
                end do

            end do

            !Skip stars
            read (fid, *)

        end do

        ! Transform THz to Hartree
        freq_ph = freq_ph/hartree_to_thz

        ! Optional output
        if (present(evec_ph_out)) evec_ph_out = evec_ph_local

        close (fid)

        ! Check if number of atoms was correct
        call terminate_if_false(mpiglobal, ias == natmtot, &
                                error_message=' Total number of atoms not &
                                correct when reading: '//trim(fname))

    end subroutine qe_read_phonon

    !> Read dielectric tensor and Born effective charges
    !> from exciting output
    subroutine exc_read_eps_inf_zstar(natoms, fname, eps_infty, zstar)

        use m_getunit, only: getunit
        use errors_warnings, only: terminate_if_false
        use modmpi, only: mpiinfo

        ! I/O
        !> number of atoms
        integer(sp), intent(in) :: natoms(:)
        !> filename
        character(len=40), intent(in) :: fname
        !> High frequency dielectric tensor \varepsilon_\infty
        real(dp), intent(out)   ::  eps_infty(:, :)
        !> Born effective charge tensor
        real(dp), intent(out)   :: zstar(:, :, :)

        ! Help variables
        type(mpiinfo) :: mpiglobal
        !> fileindex
        integer :: un
        !> status
        integer :: stat
        !> Running index combined atoms and species
        integer(sp)  :: ias
        !> check index
        integer(sp) :: ilines
        !> Number of species
        integer(sp) :: nspecies
        !> Total number of atoms
        integer(sp) :: natmtot
        !> Running index atoms
        integer(sp)  :: iatom
        !> Running index species
        integer(sp)  :: ispecies

        ! Get total number of atoms
        natmtot = size(zstar, dim=3)

        ! Get number of species
        nspecies = size(natoms)
        call getunit(un)
        open (unit=un, file=trim(fname), status='old', action='read', iostat=stat)
        if (stat /= 0) then
            write (*, *) "Error opening file, iostat:", stat, trim(fname)
            call terminate_mpi_env(mpiglobal)
        end if

        ias = 0

        ! Skip dielectric tensor line
        read (un, *)

        do ilines = 1, 3
            read (un, *) eps_infty(ilines, :)
        end do

        ! Skip effective charges and empty line
        read (un, *)

        ! Read Born effective charge tensors
        ias = 0
        do ispecies = 1, nspecies
            do iatom = 1, natoms(ispecies)
                ias = ias + 1

                ! Skip empty line andatom info
                read (un, *)
                read (un, *)

                do ilines = 1, 3
                    read (un, *) zstar(ilines, :, ias)
                end do

            end do
        end do

        close (un)

        ! Check if number of atoms was correct
        call terminate_if_false(mpiglobal, ias == natmtot, &
                                error_message=' Total number of atoms not &
                                correct when reading: '//trim(fname))

    end subroutine exc_read_eps_inf_zstar

    !> Read phonon frequencies and eigenvectors
    !> from exciting output of the type PHONON.OUT
    subroutine exc_read_phonon(natoms, fname, freq_ph, evec_ph_out)

        use m_getunit, only: getunit
        use errors_warnings, only: terminate_if_false
        use modmpi, only: mpiinfo

        ! I/O

        !> Number of atoms per species
        integer(sp), intent(in) :: natoms(:)
        !> filename
        character, intent(in) :: fname
        !> phonon frequencies
        real(dp), intent(out)  :: freq_ph(:, :)
        !> phonon eigenvectors for output
        complex(dp), intent(out), optional  :: evec_ph_out(:, :, :, :)
        !> phonon eigenvectors for output
        complex(dp), allocatable  :: evec_ph_local(:, :, :, :)

        ! Help variables
        !> Number of qpoints
        integer(sp) :: nqpoints
        !> Total number of atoms
        integer(sp) :: natmtot
        !> Number of species
        integer(sp) :: nspecies
        !> status of reading file
        integer(sp) ::   stat
        !> Running  index cartesian coord
        integer(sp) ::  i
        !> Running  index cartesian coord
        integer(sp) ::  j
        !> Running index species
        integer(sp) ::   js
        !> Running index atoms
        integer(sp) ::    ja
        !> Running index qpoint
        integer(sp) ::  iq
        !> fileindex
        integer(sp) :: fid
        !> Running index atoms
        integer(sp)  :: iatom
        !> Running index species
        integer(sp)  :: ispecies
        !> Running index combined atoms and species
        integer(sp)  :: ias
        !> Running index phonon mode
        integer(sp) ::  imode
        !> qpoint info
        real(dp) :: vqlcheck(3)
        !> mpiinfo
        type(mpiinfo) :: mpiglobal
        !> Number of phonon modes
        integer(sp) :: n_phonon_modes

        ! Get some infos
        natmtot = size(freq_ph, dim=1)/3._dp
        n_phonon_modes = 3*natmtot
        nspecies = size(natoms)
        nqpoints = size(freq_ph, dim=2)
        allocate (evec_ph_local(3, natmtot, n_phonon_modes, nqpoints))

        call getunit(fid)
        open (unit=fid, file=trim(fname), status='old', action='read', iostat=stat)
        if (stat /= 0) then
            write (*, *) "Error opening file, iostat:", stat, fname
            call terminate_mpi_env(mpiglobal)
        end if

        do iq = 1, nqpoints
            read (fid, *) i, vqlcheck

            ! skip empty line
            read (fid, *)

            do imode = 1, n_phonon_modes
                ias = 0
                ! read frequencies
                read (fid, *) i, freq_ph(imode, iq)
                do ispecies = 1, nspecies
                do iatom = 1, natoms(ispecies)
                    ias = ias + 1
                    do i = 1, 3

                        ! read eigenvectors
                        read (fid, *) js, ja, j, evec_ph_local(i, ias, imode, iq)
                    end do
                end do
                end do

                ! skip empty line
                read (fid, *)
            end do
        end do

        close (fid)

        ! Optional output
        if (present(evec_ph_out)) evec_ph_out = evec_ph_local

        ! Check if number of atoms was correct
        call terminate_if_false(mpiglobal, ias == natmtot, &
                                error_message=' Total number of atoms not &
                                correct when reading: '//trim(fname))

    end subroutine exc_read_phonon

    !> Reads Fouier coefficients of scr Coul. int. from file
    !> I.e. it reads \( W_{GG'} \) for given qpoint and phonon mode
    subroutine read_W_ph(filename, iq, qvec, imode, freq_ph, W_ph)

        use m_getunit, only: getunit
        use math_utils, only: all_close
        use modmpi, only: mpiglobal, terminate_mpi_env

        !> q-point index (requested)
        integer(sp), intent(in) :: iq
        !> q-vector in cartesian coordinates (requested)
        real(dp), intent(in) :: qvec(3)
        !> phonon mode index (requested)
        integer(sp), intent(in) :: imode
        !> Phonon frequency for given mode and q-vector (requested)
        real(dp), intent(in) :: freq_ph
        !>filename
        character(len=40), intent(in) :: filename
        !> Fourier coeffs. screened Coulomb int.
        complex(dp), intent(out) :: W_ph(:, :)

        !> record length fourier coeffs. screened Coulomb int.
        integer(dp) :: reclen
        !> read/open status
        integer(dp) :: stat
        !>file ID
        integer :: fid
        !> Number of (G+q)-vectors (requested)
        integer(sp) :: n_gqvecs

        !> q-point index (read from file)
        integer(sp) :: iq_r
        !> q-vector in cartesian coordinates (read from file)
        real(dp) :: qvec_r(3)
        !> phonon mode index (read from file)
        integer(sp) :: imode_r
        !> Phonon frequency for given mode and q-vector (read from file)
        real(dp) :: freq_ph_r
        !> Number of (G+q)-vectors (read from file)
        integer(sp) :: n_gqvecs_r

        n_gqvecs = size(W_ph, dim=1)

        !Get record length
        inquire (iolength=reclen) iq_r, qvec_r, n_gqvecs_r, &
            imode_r, freq_ph_r, W_ph

        call getunit(fid)
        open (unit=fid, file=filename, status='old', form='unformatted',&
                & action='read', access='direct', recl=reclen, iostat=stat)

        if (stat /= 0) then
            write (*, *) "Error opening file, iostat:", stat, filename
            call terminate_mpi_env(mpiglobal)
        end if

        read (fid, rec=imode, iostat=stat) iq_r, qvec_r, n_gqvecs_r, &
            imode_r, freq_ph_r, W_ph

        if (stat /= 0) then
            write (*, *) "Error reading file, iostat:", stat, filename
            call terminate_mpi_env(mpiglobal)
        end if
        close (fid)

        ! Check consistency of requested data with saved data
        if (n_gqvecs_r /= n_gqvecs &
            .or. .not. all_close(qvec_r, qvec) &
            .or. iq_r /= iq .or. imode_r /= imode &
            .or. .not. all_close(freq_ph_r, freq_ph)) then

            write (*, '(a)') 'Error(read_W_ph):&
              & differring parameters for matrix elements (current/file): '
            write (*, '(a, 2i6)') 'ngq', n_gqvecs, n_gqvecs_r
            write (*, '(a, 3f12.6, a, 3f12.6)') 'qvec', qvec, ', ', qvec_r
            write (*, '(a, 2i6)') 'for q-point :', iq, iq_r
            write (*, '(a, 2i6)') 'phonon mode :', imode, imode_r
            write (*, '(a, 2E13.6)') 'phonon freq: ', freq_ph, freq_ph_r
            write (*, '(a)') ' file: ', trim(filename)
            call terminate_mpi_env(mpiglobal)
        end if

    end subroutine read_W_ph

    !> Writes Fouier coefficients of scr Coul. int. to file
    !> I.e. it writes W_{GG'} for given qpoint and phonon mode
    subroutine write_W_ph(filename, iq, qvec, imode, freq_ph, W_ph)

        use m_getunit, only: getunit

        !> q-point index
        integer(sp), intent(in) :: iq
        !> q-vector in cartesian coordinates
        real(dp), intent(in) :: qvec(3)
        !> phonon mode index
        integer(sp), intent(in) :: imode
        !> Phonon frequency for given mode and q-vector
        real(dp), intent(in) :: freq_ph
        !> Fourier coeffs. screened Coulomb int.
        complex(dp), intent(in) :: W_ph(:, :)
        !>filename
        character(len=40), intent(in) :: filename
        !> record length fourier coeffs. screened Coulomb int.
        integer(dp) :: reclen
        !> read/open status
        integer(dp) :: stat
        !>file ID
        integer(sp) :: fid
        !> Number of (G+q)-vectors
        integer(sp) :: n_gqvecs

        n_gqvecs = size(W_ph, dim=1)

        inquire (iolength=reclen) iq, qvec, n_gqvecs, imode, freq_ph, W_ph

        call getunit(fid)
        open (unit=fid, file=filename, form='unformatted',&
                & action='write', access='direct', recl=reclen, iostat=stat)
        if (stat /= 0) then
            write (*, *) "Error opening file, iostat:", stat, filename
            call terminate_mpi_env(mpiglobal)
        end if

        write (fid, rec=imode, iostat=stat) iq, qvec, n_gqvecs, &
            imode, freq_ph, W_ph

        if (stat /= 0) then
            write (*, *) "Error writing file, iostat:", stat, filename
            call terminate_mpi_env(mpiglobal)
        end if
        close (fid)

    end subroutine write_W_ph

    !> Checks if file is present and terminates if not
    subroutine check_file_existence(filename)
        use errors_warnings, only: terminate_if_false
        character(256), intent(in) :: filename
        logical :: file_exists

        inquire (file=trim(filename), exist=file_exists)
        call terminate_if_false(mpiglobal, file_exists, error_message='The following input file is not present: '//trim(filename))

    end subroutine check_file_existence

end module phonon_screening_file_interface
