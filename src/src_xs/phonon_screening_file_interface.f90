module phonon_screening_file_interface

    use precision, only: sp, dp, str_1024
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
    !> from Quantum Espresso matdyn.x output (xxx.eig)
    !> xxxx.eig is the dynamical matrix at the Gamma point
    subroutine qe_read_phonon(alat_qe, qvecs_in, natoms, fname, freq_ph, evec_ph_out)

        use m_getunit, only: getunit
        use unit_conversion, only: hartree_to_thz
        use errors_warnings, only: terminate_if_false
        use modmpi, only: mpiinfo
        use math_utils, only: all_close
        use constants, only: twopi

        ! I/O
        !> Number of atoms per species
        integer(sp), intent(in) :: natoms(:)
        !> filename
        character(len=40), intent(in) :: fname
        !> qvecs for which to read
        real(dp), intent(in)  :: qvecs_in(:, :)
        !> Lattice constant as used by quantum espresso
        real(dp), intent(in)  :: alat_qe
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
        !> Running index combined atoms and species
        integer(sp)  :: ias
        !> fileindex
        integer(sp) :: fid
        !> Line buffers
        character(len=str_1024) :: line, fstr 
        !> I/O status
        integer(sp) :: ios
        !> helper array phonon eigenvectors
        real(dp) :: e_real(3)
        !> helper array phonon eigenvectors
        real(dp) :: e_imag(3)
        !> mpiinfo
        type(mpiinfo) :: mpiglobal
        !> Number of phonon modes
        integer(sp) :: n_phonon_modes
        !> q-vectors as read from file (for comparison with expected q-vectors)
        real(dp), allocatable :: qvecs_read(:, :)
        !> Positions in the buffers
        integer(sp) :: pos1, pos2



        ! Get some infos
        natmtot = size(freq_ph, dim=1)/3._dp
        nspecies = size(natoms)
        nqpoints = size(freq_ph, dim=2)
        n_phonon_modes = 3*natmtot

        allocate(qvecs_read(3, nqpoints))
        allocate (evec_ph_local(3, natmtot, n_phonon_modes, nqpoints))

        call getunit(fid)
        open (unit=fid, file=trim(fname), status='old', action='read')

        iq = 0
        ! Iterate over the file
        do 
            ! Try to read the line, if EOS exit
            read(fid, '(A)', iostat=ios) line
            if (ios /= 0) exit

            ! Detect q-point
            if (index(line, 'q =') /= 0) then
                iq = iq + 1
                imode = 0
                pos1 = index(line,'=') + 1
                fstr = adjustl(line(pos1:))
                read(fstr,*) qvecs_read(:, iq)
                cycle
            end if

            ! Detect frequency line
            if (index(line, 'freq') /= 0) then
                imode = imode + 1
                ias = 0
                pos1 = index(line, '=') + 1
                pos2 = index(line(pos1:), '[') - 2
                fstr = adjustl(line(pos1:pos1+pos2))
                read(fstr, *) freq_ph(imode, iq)
                cycle
            end if

            ! Detect eigenvector line
            if (index(line, '(') /= 0 .and. index(line, 'freq') == 0) then
                ias = ias + 1
                pos1 = index(line,'(') + 1
                pos2 = index(line,')') - 1
                fstr = adjustl(line(pos1:pos2))
                read(fstr,*) e_real(1), e_imag(1), e_real(2), e_imag(2), e_real(3), e_imag(3)
                evec_ph_local(:, ias, imode, iq) = cmplx(e_real, e_imag, kind=dp)
            end if

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

        ! Check if read qvecs match expected ones
        call terminate_if_false(mpiglobal, all_close(qvecs_read * twopi / alat_qe, qvecs_in, tol=1e-3_dp), &
            error_message=' q-vecs read from file differ from expected q-vecs: '//trim(fname))

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


    !> Checks if file is present and terminates if not
    subroutine check_file_existence(filename)
        use errors_warnings, only: terminate_if_false
        character(256), intent(in) :: filename
        logical :: file_exists

        inquire (file=trim(filename), exist=file_exists)
        call terminate_if_false(mpiglobal, file_exists, error_message='The following input file is not present: '//trim(filename))

    end subroutine check_file_existence

end module phonon_screening_file_interface

