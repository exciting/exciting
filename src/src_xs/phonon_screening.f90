!> Module for computation of polar phonon screening within
!> the solution of the Bethe-Salpeter-Equation. The implementation
!> (and documentation) follow Maximilian Schebek's
!> master thesis: " Impact of electron-phonon
!> interactions on excitons in polar crystals: A first-principles study".
module phonon_screening

    use precision, only: sp, dp

    implicit none

    private

    public :: phonon_screening_launcher, phonon_screening_main, &
              gen_phonon_hamiltonian, gen_phonon_screening, select_bse_energies, &
              init_phonon_hamiltonian

contains

    !> Top level routine of phonon screening in polar materials.
    !> Initialises k-, q- and (G+q)-grids.
    !> Generates output directory and checks
    !> if necessary files for phonons exist calls phonon_screening_main.
    subroutine phonon_screening_launcher
        use phonon_screening_file_interface, only: check_file_existence
        use modmpi, only: mpiglobal
        use modinput, only: input

        character(256) :: syscommand, wphdirname, zstarfname, phononfname
        logical :: file_exists

        ! Check if files for Born effective charges + diel.constant
        ! and phonons are present
        zstarfname = input%xs%phonon_screening%zstar_file
        phononfname = input%xs%phonon_screening%phonon_file
        if (mpiglobal%rank == mpiglobal%root) then
            call check_file_existence(zstarfname)
            call check_file_existence(phononfname)
        end if

        ! Generate output directory
        wphdirname = 'WPH'
        if (mpiglobal%rank == mpiglobal%root) then
            syscommand = 'test ! -e '//trim(adjustl(wphdirname))//' && mkdir '//trim(adjustl(wphdirname))
            call system(trim(adjustl(syscommand)))
        end if

        ! Initialise universal variables
        call init0
        ! Setting up k and G+k variables
        call init1
        ! Setting up q and G+q variables (task 331 therefore reduced q-points)
        call init2
        ! Initialise and compute phonon screening
        call phonon_screening_main

    end subroutine phonon_screening_launcher

    !> Computation of the phonon contribution to the direct BSE term
    !> expressed in the two-particle (electron-hole) basis: \alpha = vck.
    !> Computes for one mode with freqeuncy
    !> freq_ph at one (k,k') combination with q=k'-k.
    !> Computation for all considered electronic bands
    !>
    !> \[ H^\text{ph}_{vc\mathbf{k},v'c'\mathbf{k}'}  =
    !>    \frac{1}{V}\sum_\nu\sum_{\mathbf{GG'}}\frac{2}{\omega_{\mathbf{q}\nu}}
    !>    W_{\mathbf{q}\nu}^{\mathbf{G}\mathbf{G}'}(\mathbf{q})
    !>    M^\mathbf{G}_{cc'}(\mathbf{k,q}) [M^\mathbf{G'}_{vv'}(\mathbf{k,q})]^*\\
    !>    \times \frac{\omega_{\mathbf{q}\nu}}{2}
    !>    \Bigg[\frac{1}{\omega_{\mathbf{q}\nu}  +
    !>    \left(\varepsilon^\text{QP}_{c\mathbf{k}} -
    !>    \varepsilon^\text{QP}_{v'\mathbf{k'}}\right) - E^\lambda  } \\
    !>    +  \frac{1}{\omega_{\mathbf{q}\nu} + \left(\varepsilon^\text{QP}_{c'\mathbf{k'}}
    !>    - \varepsilon^\text{QP}_{v\mathbf{k}}\right) -  E^\lambda }\Bigg]                \]
    !>
    !> Appearing quantities:
    !> Phonon freqeuncies \[ \omega_{\mathbf{q}\nu} \]
    !> Fourier coefficients screened Coulomb interaction
    !> \[ W_{\mathbf{q}\nu}^{\mathbf{G}\mathbf{G}'}(\mathbf{q}) \]
    !>  Matrix element of the plane wave
    !> \[ M^\mathbf{G}_{mn}(\mathbf{k,q}) =
    !>    <\psi_{n\mathbf{k}}|e^{-i\mathbf{(q+G)\cdot r}}|\psi_{m\mathbf{k'}}> \]
    !> Quasi-particle energies \[ \varepsilon^\text{QP} \]
    !> Indices of conduction and valence bands \[ c,v \]
    subroutine gen_phonon_hamiltonian(one_over_vol, wfc_ph, M_planewave_uu, &
                  & M_planewave_oo, eigvals_o, eigvals_u, &
                  & freq_ph, eigen_exc, map_alpha_vck, iaoff, jaoff, H_ph)
        use constants, only: zone, zzero
        use xlapack, only: matrix_multiply

        !>  Second-variational eigenvalues of considered occupied bands at k,k'
        real(dp), intent(in) :: eigvals_o(:, :)
        !>  Second-variational eigenvalues of considered unoccupied bands at k,k'
        real(dp), intent(in) :: eigvals_u(:, :)
        !> 1/Vol, Vol = BvK Volume: nkpt*V_{unit_cell}
        complex(8), intent(in) ::  one_over_vol
        !> Plane-wave matrix elements between occupied (= valence) states
        complex(dp), intent(in) ::  M_planewave_oo(:, :)
        !> Plane-wave matrix elements between unoccupied (= conduction) states
        complex(dp), intent(in) ::  M_planewave_uu(:, :)
        !> Fourier coefficients of phonon screended Coulomb interation W_{GG'}
        !> for one q_point and one phonon mode scaled by 2/freq_ph
        complex(8), intent(in) :: wfc_ph(:, :)
        !> Phonon frequency of considered mode at q=k'-k
        real(dp), intent(in) :: freq_ph
        !> Excitation energy of considered excitation
        real(dp), intent(in) :: eigen_exc
        !> Mapping between alpha and relative vck
        integer(sp), intent(in) :: map_alpha_vck(:, :)
        !> Offset for bands
        integer(sp), intent(in):: iaoff, jaoff
        !> ik,jk block of phonon direct BSE hamiltonian in two-particle basis
        complex(dp), intent(out):: H_ph(:, :)

        ! Help arrays
        !> Number of unoccupied bands squared
        integer(sp) :: n_unocc_bands_sqrd
        !> Number of occupied bands squared
        integer(sp) :: n_occ_bands_sqrd
        !> Number of (G+q)-vectors for current (k,k')
        integer(sp) :: n_gq_vectors
        !> Intermediate array for triple matrix multiplication
        complex(dp), allocatable :: zm(:, :)
        !> Intermediate BSE term with indices (vv', cc')
        complex(dp), allocatable :: H_ph_comb(:, :)
        !> Elements of the  BSE term resolved into individual bands (v,c, v', c')
        complex(8), allocatable :: H_ph_bands(:, :, :, :)
        !> Number of participating occupied bands
        integer(sp) :: n_occ_bands
        !> Number of participating unoccupied bands
        integer(sp) :: n_unocc_bands
        !> Runing index occupied states
        integer(sp) :: io, jo
        !> Runing index unoccupied states
        integer(sp) :: iu, ju
        !> Combined indices
        integer(sp) :: j1, j2, ia, ja
        !> Energy differences for bands (ck - v'k')  and (c'k' - vk)
        real(dp) :: energy_diff_1, energy_diff_2
        !> weight of mode
        real(dp) :: mode_weight
        !> Eigenenergies of states at k (ki) and k' (kj)
        real(dp) :: eigen_c_ki, eigen_c_kj, eigen_v_ki, eigen_v_kj

        ! Get sizes and allocate arrays
        n_gq_vectors = size(M_planewave_oo, dim=2)
        n_occ_bands_sqrd = size(M_planewave_oo, dim=1)
        n_unocc_bands_sqrd = size(M_planewave_uu, dim=1)
        n_occ_bands = int(n_occ_bands_sqrd**(0.5_dp))
        n_unocc_bands = int(n_unocc_bands_sqrd**(0.5_dp))
        allocate (zm(n_unocc_bands_sqrd, n_gq_vectors))
        allocate (H_ph_comb(n_unocc_bands_sqrd, n_occ_bands_sqrd))
        allocate (H_ph_bands(n_unocc_bands, n_occ_bands, n_unocc_bands, n_occ_bands))

        ! H_ph_comb = 1/V \sum_{GG'} M^{G',*}_{vv'} M^G_{cc'} W^{ph}_{GG'}
        call matrix_multiply(M_planewave_uu, wfc_ph, zm, 'n', 'n')
        call matrix_multiply(zm, M_planewave_oo, H_ph_comb, 'n', 't')
        H_ph_comb = H_ph_comb*one_over_vol

        ! The current matrix H_ph_comb has indices (vv', cc'). To get to indices
        ! (vc, v'c') we need to resolve into individual band contributions
        ! and build the combined index afterwards.

        ! Map back to individual band indices
        do jo = 1, n_occ_bands
            eigen_v_kj = eigvals_o(jo, 2)
            do ju = 1, n_unocc_bands
                eigen_c_kj = eigvals_u(ju, 2)
                do io = 1, n_occ_bands
                    eigen_v_ki = eigvals_o(io, 1)
                    do iu = 1, n_unocc_bands
                        eigen_c_ki = eigvals_u(iu, 1)

                        j1 = iu + (ju - 1)*n_unocc_bands
                        j2 = io + (jo - 1)*n_occ_bands
                        ! Ht_{iu_j1 ju_j1, io_j2 jo_j2} -> H_{iu_j1 io_j2, ju_j1 jo_j2}
                        H_ph_bands(iu, io, ju, jo) = H_ph_comb(j1, j2)

                        !Weight each element
                        energy_diff_1 = eigen_c_ki - eigen_v_kj
                        energy_diff_2 = eigen_c_kj - eigen_v_ki
                        mode_weight = -freq_ph/2._dp* &
                                      (1._dp/(freq_ph + energy_diff_1 - eigen_exc) + &
                                       1._dp/(freq_ph + energy_diff_2 - eigen_exc))

                        H_ph_bands(iu, io, ju, jo) = mode_weight*H_ph_bands(iu, io, ju, jo)
                    end do
                end do
            end do
        end do

        ! Compute mode contribution to BSE Hamiltonian in two-particle basis.
        ! The composite index corresponds to one independent-particle
        ! transition, i.e. H_ph has indices (vck, v'c'k')
        do ja = 1, n_occ_bands*n_unocc_bands
            do ia = 1, n_occ_bands*n_unocc_bands
                ju = map_alpha_vck(1, ja + jaoff)
                jo = map_alpha_vck(2, ja + jaoff)
                iu = map_alpha_vck(1, ia + iaoff)
                io = map_alpha_vck(2, ia + iaoff)
                ! Ht_{iu_j1 io_j2, ju_j1 jo_j2} -> Hab_{iu_a io_a, ju_b jo_b}
                H_ph(ia, ja) = H_ph_bands(iu, io, ju, jo)
            end do
        end do

    end subroutine gen_phonon_hamiltonian

    !> Reads from file the phonon screened Coulomb interaction W_{GG'}(q,nu)
    !> and the phonon frequencies for all q-vectors and all modes.
    subroutine init_phonon_hamiltonian(n_phonon_modes, n_atoms_spec, n_gq, &
                                       qvecs, scieffg_ph, freq_ph)
        use modinput, only: input
        use constants, only: zzero
        use phonon_screening_file_interface, only: qe_read_phonon, &
                                           read_W_ph, exc_read_phonon

        !> Number of phonon modes
        integer(sp), intent(in) :: n_phonon_modes
        !> Number atoms per species
        integer(sp), intent(in) :: n_atoms_spec(:)
        !> Number of (G+q)-vectors for reduced q-vectors
        integer(sp), intent(in) :: n_gq(:)
        !> Reduced set of q-vectors
        real(dp), intent(in) :: qvecs(:, :)

        !> Phonon frequencies for all modes and all q-vectors
        real(dp), intent(out) :: freq_ph(:, :)
        !> Fourier coefficients of the phonon screened Coulomb W_{GG'}(q,nu)
        !> for all q-vectors for and all phonon modes
        complex(dp), intent(out) :: scieffg_ph(:, :, :, :)

        !> Number of reduced q-vectors
        integer(sp) :: nqptr
        !> filename
        character(len=40) :: fname
        !> String to number
        character(256) :: strnum
        !> Running index reduced q-vectors
        integer(sp) :: iqr
        !> Running index non-reduced q-vectors
        integer(sp) :: iqrnr
        !> Running index phonon modes
        integer(sp) :: imode

        scieffg_ph = zzero

        nqptr = size(n_gq)

        ! Read phonon frequencies. Check if file comes from QE or exciting
        fname = input%xs%phonon_screening%phonon_file
        if (input%xs%phonon_screening%file_type == "quantum_espresso") then
            call qe_read_phonon(freq_ph=freq_ph, natoms=n_atoms_spec, fname=fname)
        elseif (input%xs%phonon_screening%file_type == "exciting") then
            call exc_read_phonon(freq_ph=freq_ph, natoms=n_atoms_spec, fname=fname)
        end if

        do iqr = 1, nqptr

            ! Read phonon screened Coulomb interaction
            write (strnum, '(I5.5)') iqr
            fname = "WPH/wph_Q"//trim(strnum)//".OUT"
            do imode = 1, n_phonon_modes
                call read_W_ph(fname, iqr, qvecs(:, iqr), imode, &
                               freq_ph(imode, iqr), &
                               scieffg_ph(1:n_gq(iqr), 1:n_gq(iqr), imode, iqr))
            end do

        end do

    end subroutine

    !> Reads Quantum Espresso or exciting data for Born charges, phonon
    !> frequs and diel. tensor from file and stores them in a consistent format.
    !> Computes and writes to file the phonon screened Coulomb interaction
    !> \[ W_{\mathbf{q}\nu}^{\mathbf{G}\mathbf{G}'}(\mathbf{q}) \]
    !> for all q-points and all phonon modes.
    subroutine phonon_screening_main

        use constants, only: zzero
        use modxs, only: ngq, qpari, qparf, vgqc
        use mod_qpoint, only: nqpt, vqc
        use mod_atoms, only: natmtot, natoms, spmass, atposc
        use mod_lattice, only: omega
        use mod_kpoint, only: nkpt
        use modinput, only: input
        use exciting_mpi, only: xmpi_bcast
        use modmpi, only: mpiglobal
        use phonon_screening_file_interface, only: qe_read_eps_inf_zstar, qe_read_phonon, &
        & write_W_ph, exc_read_eps_inf_zstar, exc_read_phonon
        use math_utils, only: all_zero
        !> phonon frequencies
        real(dp), allocatable    :: freq_ph(:, :)
        !> phonon eigenvectors
        complex(dp), allocatable  :: evec_ph(:, :, :, :)
        !> High frequency dielectric tensor \[ \varepsilon_\infty \]
        real(dp), allocatable:: eps_infty(:, :)
        !> Born effective charge tensor
        real(dp), allocatable:: zstar(:, :, :)
        !> Running index qpoint
        integer :: iq
        !> Running index phonon modes
        integer ::  imode
        !> Number of phonon modes
        integer(sp) :: n_phonon_modes
        !> Vector of direction for approaching q=0
        real(dp) :: qdir(3)
        !> filenames
        character(len=40) :: phonon_fname, zstar_fname, wph_fname
        !> string for filename
        character(256) :: strnum
        !> Tolerance frequency for neglecting phonon mode
        real(dp), parameter :: freq_tol = 1.E-5_dp
        !> Fourier coefficients of phonon screended Coulomb interation W_{GG'}
        !> for one q_point and one phonon mode
        complex(dp), allocatable   ::  W_ph(:, :)

        n_phonon_modes = 3*natmtot

        !Allocate arrays on all procs and set zero
        allocate (freq_ph(n_phonon_modes, nqpt))
        allocate (evec_ph(3, natmtot, n_phonon_modes, nqpt))
        allocate (eps_infty(3, 3))
        allocate (zstar(3, 3, natmtot))

        if (mpiglobal%rank == mpiglobal%root) then
            zstar_fname = input%xs%phonon_screening%zstar_file
            phonon_fname = input%xs%phonon_screening%phonon_file

            ! Check if files come from exciting or QE calculation
            if (input%xs%phonon_screening%file_type == "quantum_espresso") then
                call qe_read_eps_inf_zstar(natoms, zstar_fname, &
                                           eps_infty, zstar)
                call qe_read_phonon(natoms, phonon_fname, freq_ph, evec_ph)
            elseif (input%xs%phonon_screening%file_type == "exciting") then
                call exc_read_eps_inf_zstar(natoms, zstar_fname, &
                                            eps_infty, zstar)
                call exc_read_phonon(natoms, phonon_fname, freq_ph, evec_ph)
            end if

        end if

        Call xmpi_bcast(mpiglobal, eps_infty)
        Call xmpi_bcast(mpiglobal, zstar)
        Call xmpi_bcast(mpiglobal, freq_ph)
        Call xmpi_bcast(mpiglobal, evec_ph)

        ! Direction in reciprocal space which was used to
        ! evaluate the polar phonons at Gamma. In Quantum Espresso, this is
        ! the direction from the second q-vector in the list to Gamma
        qdir = vqc(1:3, 2)/norm2(vqc(:, 2))

        ! Distribute reduced q-points
        call genparidxran('q', nqpt)

        do iq = qpari, qparf

            !Allocate final quantity W_{GG'} for given phonon mode and qpoint
            if (allocated(W_ph)) deallocate (W_ph)
            allocate (W_ph(ngq(iq), ngq(iq)))

            !Output file
            write (strnum, "(I5.5)") iq
            wph_fname = "WPH/wph_Q"//trim(strnum)//".OUT"

            do imode = 1, n_phonon_modes
                W_ph = zzero

                ! If phonon frequency is too small just put zero
                if (.not. all_zero(freq_ph(imode, iq), freq_tol)) then

                    !Compute final quantity W_{GG'} for given phonon mode and qpoint
                    call gen_phonon_screening(qdir, vgqc(:, :ngq(iq), iq), eps_infty, zstar, &
                                    & freq_ph(imode, iq), evec_ph(:, :, imode, iq), nkpt, omega, natoms, spmass, atposc, W_ph)

                end if

                call write_W_ph(wph_fname, iq, vqc(:, iq), imode, &
                                freq_ph(imode, iq), W_ph)
            end do
        end do

    end subroutine phonon_screening_main

    !> Computes Fourier coeff. screened Coulomb int. due to phonons in polar materials
    !> For non-polar materials: All coefficients vanish due to
    !> zero Born-effective charges.
    !> Computation for one phonon mode and and one q-point (Eq. 4.13b):
    !>
    !> \[ W_{\mathbf{q}\nu}^{\mathbf{G}\mathbf{G}'}(\mathbf{q}) =
    !>   \tilde{g}_{\mathbf{q}\nu}^{\mathbf{G}}(\mathbf{q})
    !>   \left[\tilde{g}_{\mathbf{q}\nu}^{\mathbf{G}'}(\mathbf{q})\right]^* \]
    !>
    !> where \[ \tilde{g}_{\mathbf{q}\nu}^{\mathbf{G}}(\mathbf{q}) \] is the
    !> Fourier coefficient of the electron-phonon coupling function.
    !> TODO (Max) Issue 59: Collect data in types -> pass less arguments
    subroutine gen_phonon_screening(qdir, gqvecs, eps_infty, zstar, &
      & freq_ph, evec_ph, nkpt, omega, n_atom_species, spmass, atposc, W_ph)

        use constants, only: zzero, fourpi, twopi, pi, zi, zone
        use math_utils, only: all_zero
        use xlapack, only: matrix_multiply, outer_product
        use asserts, only: assert
        use modmpi, only: mpiglobal, terminate_mpi_env

        !> High frequency dielectric tensor \varepsilon_\infty
        real(dp), intent(in):: eps_infty(:, :)
        !> Born effective charge tensor
        real(dp), intent(in):: zstar(:, :, :)
        !> frequency of phonon mode at q-point
        real(dp), intent(in)    ::  freq_ph
        !> phonon eigenvectors of phonon mode at q-point
        complex(dp), intent(in)  :: evec_ph(:, :)
        !> Unit cell volume
        real(dp), intent(in) :: omega
        !> (q+G)-vectors for q-point
        real(dp), intent(in) :: gqvecs(:, :)
        !> Unit vector for approaching Gamma
        real(dp), intent(in) :: qdir(3)
        !> number of atoms per species
        integer(sp), intent(in) ::  n_atom_species(:)
        !> number of k-points (different from number of reduced q-points)
        integer(sp), intent(in) :: nkpt
        !> Species masses
        real(dp), intent(in):: spmass(:)
        !> Atomic positions in cartesian coordinates
        Real(dp), intent(in) :: atposc(:, :, :)
        !> Fourier coefficients of phonon screended Coulomb interation W_{GG'}
        !> for one q_point and one phonon mode
        complex(dp), intent(out)   ::  W_ph(:, :)
        !> Indeces of finite G+q
        integer(sp), allocatable :: indices_finite_gq(:)
        !> Index of G+q = 0
        integer(sp), allocatable :: indices_zero_gq(:)

        ! Help variables
        !> Fourier coefficent coupling function, nedded to construct W_ph
        complex(dp), allocatable :: coupl_fourier(:)

        !> Polarization vector of a mode, weighted by exp atomic positions
        complex(dp), allocatable :: pol_vec(:, :)
        !> number of G-vectors for q-point
        integer :: n_gqvecs
        !> Running index (q+G)point
        integer :: igq
        !> Prefactor
        complex(dp) :: prefac
        !> Running index
        integer(sp) :: i

        ! Get number of (G+q)-vectors for current q-point
        n_gqvecs = size(W_ph, dim=2)

        call assert(n_gqvecs == size(gqvecs, dim=2), &
                    "Error gen_phonon_screening: Number of (G+q)-vectors not &
                    matching the dimension screened Coulomb interaction")

        ! Find indices of zero and finite (G+q)-vectors, respectively
        indices_finite_gq = pack([(i, i=1, n_gqvecs)], norm2(gqvecs, dim=1) > 1.e-9_dp)
        indices_zero_gq = pack([(i, i=1, n_gqvecs)], norm2(gqvecs, dim=1) < 1.e-9_dp)
        if (size(indices_zero_gq) > 1) then
            call terminate_mpi_env(mpiglobal, &
                                   'Error gen_phonon_screening: More than one (G+q) &
                                   -vector is zero, but should be only the case for G=q=0')
        end if

        allocate (coupl_fourier(n_gqvecs))
        allocate (pol_vec(3, n_gqvecs))

        prefac = fourpi/omega*zi/sqrt(2._dp*abs(freq_ph))

        call calc_weighted_polarization_vec(zstar, spmass, atposc, &
                                            evec_ph, gqvecs(:, 1:n_gqvecs), &
                                            n_atom_species, pol_vec)

        do i = 1, size(indices_finite_gq)
            igq = indices_finite_gq(i)
            call calc_coupling_fourier(prefac, gqvecs(:, igq), &
                                       pol_vec(:, igq), eps_infty, coupl_fourier(igq))
        end do

        do i = 1, size(indices_zero_gq)
            igq = indices_zero_gq(i)
            ! qdir unit vector for approaching q=0
            call calc_coupling_fourier(prefac, qdir, pol_vec(:, igq), eps_infty, coupl_fourier(igq))
        end do

        ! Build 2/freq_ph*W_{GG'} for given q-point and phonon mode
        ! Extension by 2/freq_ph corresponds to static screening
        ! Inclusion of dynamics using a weight in (0,1) when constructing phonon
        ! BSE-Hamiltonian in gen_phonon_hamiltonian
        W_ph = zzero
        call outer_product(coupl_fourier, coupl_fourier, W_ph, conjg_b=.true.)
        W_ph = W_ph*2._dp/freq_ph*omega

        ! Analytically deal with divergence of head and wings for q=0
        if (size(indices_zero_gq) > 0) then
            call average_head_wings(indices_finite_gq, indices_zero_gq, W_ph, nkpt*omega)
        end if

    end subroutine gen_phonon_screening

    !> Analytically deals with the divergences of  head and wings of
    !> the phonon screened Coulomb interaction for q=0 by means of
    !> a spherical average. The head decays as 1/q^2,  wings as 1/q,
    !> and the corresponding elements are averaged over a subcell in
    !> reciprocal space with radius \[ q_s = \left(6\pi^2/V\right)^{1/3} ]
    !> and volume \[  V_q = 2\pi^3/V \], where V is the
    !> crystal volume = N_kpt * unit cell volume.
    !> The results for the diverging factors are (Eq. 4.43 - 44):
    !>
    !> \[ \frac{1}{V_q}\int_{V_q}\text{d}\mathbf{p}\,
    !>    \frac{1}{|\mathbf{p}|^2} =  \frac{4\pi}{V_q}\,q_s ],
    !>
    !> \[ \frac{1}{V_q}\int_{V_q}\text{d}\mathbf{p}\,\frac{1}{|\mathbf{p}|}
    !>    = \frac{2\pi}{V_q}\,q^2_s \]
    subroutine average_head_wings(ids_fin, ids_zero, W_ph, vol)
        use constants, only: pi, twopi, fourpi
        integer(sp), intent(in) :: ids_fin(:)
        integer(sp), intent(in) :: ids_zero(:)
        !> Fourier coefficients of phonon screended Coulomb interation W_{GG'}
        !> for one q_point and one phonon mode
        complex(dp), intent(inout)   ::  W_ph(:, :)
        !> Volume of sphere in rec. space for integrating out divergence
        real(dp) :: volq
        !> Radius of sphere in rec. space for integrating out divergence
        real(dp) :: radius_q_singularity
        !> Born von Karman Volume nkpt*omega (omega is unit cell volume)
        real(dp), intent(in) :: vol
        !> Index of (G+q)=0 vector
        integer(sp) :: id

        id = ids_zero(1)
        volq = (twopi**3)/vol
        radius_q_singularity = (6*(pi**2)/vol)**(1._dp/3._dp)

        W_ph(id, id) = W_ph(id, id)*(1/volq)*radius_q_singularity*fourpi
        W_ph(id, ids_fin) = W_ph(id, ids_fin)*(1/volq)*twopi*radius_q_singularity**2
        W_ph(ids_fin, id) = W_ph(ids_fin, id)*(1/volq)*twopi*radius_q_singularity**2

    end subroutine

    !> This routine computes the Fourier coefficient of the
    !> electron-phonon coupling function for one (G+q)-vector,
    !> and one phonon frequency. The coefficient is defined by (Eq. 4.33)
    !>
    !> \[  \Tilde{g}_{\mathbf{q}\nu}^{\mathbf{G}}(\mathbf{q}) =
    !>     i\, B_{\mathbf{q}\nu}\,\mathbf{C(q+G)}
    !>     \cdot\bar{\mathbf{Q}}_{\nu\mathbf{G}}(\mathbf{q})    \],
    !> where the following definitions were applied (Eqs. 4.28a-d):
    !>
    !> \[ B_{\mathbf{q}\nu}  = i \frac{4\pi}{\Omega_\text{UC}}
    !>                  \left( \frac{1}{2\omega_{\mathbf{q}\nu}} \],
    !>
    !> \[ \mathbf{C(q+G)}  = \frac{\mathbf{(q+G)}}{\mathbf{(q+G)}
    !>    \cdot\boldsymbol{\varepsilon}_\infty\cdot\mathbf{(q+G)}}
    !>    = \frac{\mathbf{X}_\mathbf{G}
    !>        (\boldsymbol{\varepsilon}_\infty,\mathbf{q})}
    !>         {|\mathbf{q+G}|}                                             \],
    !>
    !> \[ \mathbf{X}_\mathbf{G}(\boldsymbol{\varepsilon}_\infty,\mathbf{q})
    !>    = \frac{\hat{e}_\mathbf{q+G}}{\hat{e}_\mathbf{q+G}\cdot
    !>      \boldsymbol{\varepsilon}_\infty\cdot\hat{e}_\mathbf{q+G}} \],
    !>
    !> \[ \bar{\mathbf{Q}}_{\nu\mathbf{G}}(\mathbf{q})  =
    !>    \sum_\kappa\frac{1}{\sqrt{M_\kappa}}\,\mathbf{Z}^*_\kappa
    !>    \cdot\mathbf{e}_{\kappa\nu}(\mathbf{q})\,e^{-i(\mathbf{q+G})
    !>    \cdot\boldsymbol{\tau}_{\kappa}}                            \],
    !>
    !> In the code we use
    !> \[ e^{-i(\mathbf{q+G})\cdot\boldsymbol{\tau}_{\kappa}}
    !>    =\cos[(\mathbf{q+G}) \cdot \boldsymbol{\tau}_{\kappa}]
    !>    - i\sin[(\mathbf{q+G}) \cdot \boldsymbol{\tau}_{\kappa}]    \].
    !>
    !> We used:
    !> Born-effective charge tensor \[ \mathbf{Z}^*_\kappa \]
    !> Atomic masses \[ M_\kappa \]
    !> and positions \[ \boldsymbol{\tau}_{\kappa} \]
    !> Phonon frequencies \[ \omega_{\mathbf{q}\nu} \]
    !> and eigenvectors  \[ \mathbf{e}_{\kappa\nu}(\mathbf{q})  \]
    !> Unit cell volume \[ \Omega_\text{UC}  \]
    !> High-frequency dielectric tensor \[ \boldsymbol{\varepsilon}_\infty \]
    subroutine calc_coupling_fourier(prefac, gqvec, pol_vec, &
                                     eps_infty, coupl_fourier)
        use xlapack, only: matrix_multiply
        !> Prefactor \[ B_{\mathbf{q}\nu} \]
        complex(dp), intent(in) :: prefac
        !> (G+q)-vector.
        real(dp), intent(in) :: gqvec(3)
        !> Polarization vector \[\bar{\mathbf{Q}}_{\nu\mathbf{G}}(\mathbf{q})\]
        complex(dp), intent(in) :: pol_vec(3)
        !> Fourier coefficient of the coupling function
        !> \[ \Tilde{g}_{\mathbf{q}\nu}^{\mathbf{G}}(\mathbf{q}) \]
        complex(dp), intent(out) :: coupl_fourier
        !> High frequency dielectric tensor \[ \varepsilon_\infty \]
        real(dp), intent(in):: eps_infty(:, :)

        !> Running index
        integer(sp) :: i
        !> Index of (G+q)-vector
        integer(sp) :: igq
        !> Help vector defined by matrix-vector multiplication eps_inf*(q+G)
        real(dp) :: epsinf_egq(3)
        !> matrix multiplication (q+G)\cdot\varepsilon_\infty\cdot(q+G)
        real(dp) ::  egq_epsinf_egq
        !> Unit vector in (G+q)-direction
        real(dp) :: unit_gq(3)
        !> Magnitude of (G+q)-vector
        real(dp) :: mag_gq

        mag_gq = norm2(gqvec)
        unit_gq = gqvec/mag_gq

        call matrix_multiply(eps_infty, unit_gq, epsinf_egq)
        egq_epsinf_egq = dot_product(unit_gq, epsinf_egq)
        coupl_fourier = dot_product(pol_vec, unit_gq)/egq_epsinf_egq/mag_gq
        coupl_fourier = coupl_fourier*prefac

    end subroutine

    !> Computes the weighted polarization vector for all (G+q)-vectors
    !> and one phonon mode \[ \bar{\mathbf{Q}}_{\nu\mathbf{G}}(\mathbf{q}) \]
    ! (see subroutine calc_coupling_fourier for definition).
    subroutine calc_weighted_polarization_vec(zstar, spmass, atposc, evec_ph, gqvecs, n_atom_species, pol_vec)

        use xlapack, only: matrix_multiply
        use constants, only: zzero
        !> Born effective charge tensor
        real(dp), intent(in):: zstar(:, :, :)
        !> Species masses
        real(dp), intent(in):: spmass(:)
        !> Atomic positions in cartesian coordinates
        Real(dp), intent(in) :: atposc(:, :, :)
        !> phonon eigenvectors of phonon mode at q-point
        complex(dp), intent(in)  :: evec_ph(:, :)
        !> (q+G)-vectors for q-point
        real(dp), intent(in) :: gqvecs(:, :)
        !> number of atoms per species
        integer(sp), intent(in) :: n_atom_species(:)
        !> Polarization vector of a mode, weighted by exp atomic positions
        complex(dp), intent(out) :: pol_vec(:, :)

        !> number of species
        integer(sp) :: nspecies
        !> Running index species
        integer :: ispecies
        !> Running index combined atoms and species
        integer :: ias

        !> Running index (q+G)point
        integer(sp) :: igq
        !> number of G-vectors for q-point
        integer(sp) :: n_gqvecs
        !> Born effective charge tensor times phonon eigenvectors
        complex(dp), allocatable :: zstar_evec_ph(:, :)
        !>  prefactor
        real(dp) :: t1
        !>  prefactor
        real(dp) :: t2
        !> Running index atoms
        integer :: iatom
        !> total number of atoms
        integer(sp) :: n_atom_total

        n_gqvecs = size(gqvecs, dim=2)

        ! Get number of species
        nspecies = size(n_atom_species)

        ! Get total number of atoms
        n_atom_total = size(zstar, dim=3)

        allocate (zstar_evec_ph(3, n_atom_total))

        pol_vec = zzero

        do igq = 1, n_gqvecs
            ias = 0
            do ispecies = 1, nspecies
                do iatom = 1, n_atom_species(ispecies)
                    ias = ias + 1

                    call matrix_multiply(zstar(:, :, ias), evec_ph(:, ias), zstar_evec_ph(:, ias))

                    t2 = -dot_product(gqvecs(:, igq), atposc(:, iatom, ispecies))

                    t1 = 1._dp/sqrt(spmass(ispecies))
                    pol_vec(:, igq) = pol_vec(:, igq) + &
                                      t1*cmplx(cos(t2), sin(t2), dp)*zstar_evec_ph(:, ias)
                end do
            end do

        end do

    end subroutine

    !> Takes all computed Kohn-Sham energy eigenvalues
    !> and returns those of the states which are
    !> considered in BSE calculation.
    !> Adds a scissor shift to the conduction bands.
    subroutine select_bse_energies(limits_bse, eigvals_k, eigvals_kp, scissor,&
                                    & eigvals_o_bse, eigvals_u_bse)

        !> Upper/lower bounds for valence and conduction bands
        integer(sp), intent(in) :: limits_bse(4)
        !> Scissor paratemer
        real(dp), intent(in) :: scissor
        !> All computed energy eigenvalues at k
        real(dp), intent(in) :: eigvals_k(:)
        !> All computed KS energy eigenvalues at k'
        real(dp), intent(in) :: eigvals_kp(:)
        !> Energy eigenvalues of occupied the states considered in BSE at k,k'
        real(dp), intent(out) :: eigvals_o_bse(:, :)
        !> Energy eigenvalues of unoccupied states considered in BSE at k,k'
        real(dp), intent(out) :: eigvals_u_bse(:, :)

        !> Lowest occupied state (absolut index)
        integer(sp) ::  lowest_occ_state
        !> Highest occupied state (absolut index)
        integer(sp) ::  highest_occ_state
        !> Lowest occupied state (absolut index)
        integer(sp) ::  lowest_unocc_state
        !> Lowest occupied state (absolut index)
        integer(sp) ::  highest_unocc_state

        lowest_unocc_state = limits_bse(1)
        highest_unocc_state = limits_bse(2)
        lowest_occ_state = limits_bse(3)
        highest_occ_state = limits_bse(4)

        ! Eigenvalues of occ bands at k
        eigvals_o_bse(:, 1) = eigvals_k(lowest_occ_state:highest_occ_state)
        ! Eigenvalues of occ bands at k'
        eigvals_o_bse(:, 2) = eigvals_kp(lowest_occ_state:highest_occ_state)
        ! Eigenvalues of unocc bands at k
        eigvals_u_bse(:, 1) = eigvals_k(lowest_unocc_state:highest_unocc_state) + &
                          & scissor
        ! Eigenvalues of unocc bands at k'
        eigvals_u_bse(:, 2) = eigvals_kp(lowest_unocc_state:highest_unocc_state) + &
                          & scissor

    end subroutine select_bse_energies

end module phonon_screening
