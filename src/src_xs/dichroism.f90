!> This modules contains all routines for the calculation of
!> valley-selective circular dichroism (VSCD) in 2d hexagonal materials
!> as presented in Caruso et al. J. Phys. Chem. Lett. 2022, 13, 5894−5899
!> (in the following only referred to as "Caruso et al." ).
!> The second important publication relevant to this implementation is
!> the Vorwerk et al 2019 Electron. Struct. 1 037001 (in the following only referred to as "Vorwerk et al." )
!> where the BSE implementation of exciting is documented.
module dichroism
    use precision, only: dp
    use modmpi, only: mpiglobal
    use constants, only: zi
    implicit none

    private

    public :: dichroism_main

contains

    !> Main routine, gathers all relevant data from global variables
    !> and passes it to the computing routines.
    subroutine dichroism_main(input, bse_gap, exc_eigvals, exc_eigvecs)

        use modbse, only: sci, smap, smap_rel, ensortidx, &
                          nk_bse, de, eval0, koulims
        use modxs, only: unitout, vkl0
        use mod_lattice, only: bvec, omega
        use mod_kpoint, only: vkcnr
        use modinput, only: input_type
        use mod_eigenvalue_occupancy, only: evalsv
        use m_genwgrid, only: genwgrid
        use os_utils, only: make_directory_command
        use modmpi, only: mpiglobal 
        use constants, only: pi

        !> BSE eigenvalues
        real(dp), intent(in) :: exc_eigvals(:)
        !> BSE band gap (smallest gap over all participating k-points)
        real(dp), intent(in) :: bse_gap
        !> BSE eigenvectors
        complex(dp), intent(in) :: exc_eigvecs(:, :)
        !> Input to get some global information about calculation
        type(input_type), intent(in) :: input

        !> K valley-resolved dichroic tensor
        complex(dp), allocatable :: dichroic_k(:, :, :)
        !> K_bar valley-resolved dichroic tensor
        complex(dp), allocatable :: dichroic_kbar(:, :, :)
        !> K-K_bar valley-resolved dichroic tensor
        complex(dp), allocatable:: dichroic_k_kbar(:, :, :)
        !> Number of eigenvalues found when diagonalizing BSE Hamiltonian
        integer :: n_exc
        !> Frequencies used for spectrum generation
        real(dp), allocatable :: frequencies(:)
        !> Number of frequencies used for spectrum generation
        integer :: n_freqs
        !> Size of the BSE Hamiltonian
        integer :: ham_size
        !> Momentum transfer index (always q=0)
        integer, parameter ::  iqmt = 1
        !> K valley-resolved transition coefficients
        complex(dp), allocatable :: trans_coeffs_k(:, :)
        !> K_bar valley-resolved transition coefficients strength
        complex(dp), allocatable :: trans_coeffs_kbar(:, :)
        !> Momentum matrix elements between unoccupied and occupied states
        complex(dp), allocatable :: pmat_unocc_occ(:, :, :, :)
        !> K valley-resolved dipole matrix elements
        complex(dp), allocatable :: dipole_mat_k(:, :)
        !> K_bar valley-resolved dipole matrix elements
        complex(dp), allocatable :: dipole_mat_kbar(:, :)
        !> K valley-resolved dichroic tensor
        complex(dp), allocatable :: dichroic_tensor_k(:, :, :)
        !> K_bar valley-resolved dichroic tensor
        complex(dp), allocatable ::  dichroic_tensor_kbar(:, :, :)
        !> Mixed valley-resolved dichroic tensor
        complex(dp), allocatable ::  dichroic_tensor_kkbar(:, :, :)
        !> K valley-resolved oscillator strength
        complex(dp), allocatable ::  osci_strengths_k(:, :, :)
        !> K_bar valley-resolved oscillator strength
        complex(dp), allocatable ::  osci_strengths_kbar(:, :, :)
        !> Mixed valley-resolved oscillator strength
        complex(dp), allocatable ::             osci_strengths_k_kbar(:, :, :)
        !> True if independent-particle approximation was employed in BSE
        logical :: use_ip
        !> True if BSE is run on top of GW
        logical :: use_gw
        !> Names of directories where to store dichroic tensors
        character(256), parameter :: dirnames_dichroic_tensors(3) = [character(len=256) :: &
                                                                     'DICHROIC_K', 'DICHROIC_KBAR', 'DICHROIC_K_KBAR']
        !> Names of directories where to store oscillator strengths
        character(256), parameter :: dirnames_oscillator_strengths(3) = [character(len=256) :: &
                                                                         'OSCI_K', 'OSCI_KBAR', 'OSCI_K_KBAR']
        !> Running index directories
        integer :: idir
        !> OS command
        character(256) :: os_command

        do idir = 1, size(dirnames_dichroic_tensors)
            os_command = make_directory_command(dirnames_dichroic_tensors(idir))
            call system(trim(adjustl(os_command)))
            os_command = make_directory_command(dirnames_oscillator_strengths(idir))
            call system(trim(adjustl(os_command)))
        end do

        use_gw = associated(input%gw)
        use_ip = input%xs%bse%bsetype == "IP"

        n_exc = size(exc_eigvals)

        ham_size = size(exc_eigvals)

        ! Get Momentum matrix elements from file
        call read_momentum_matrix_elements(koulims, use_gw, vkl0, evalsv, eval0, pmat_unocc_occ)

        ! Build valley-resolved dipole matrix elements
        call calc_k_resolved_dipole_matrix_elements(pmat_unocc_occ, smap, smap_rel, &
                                   vkcnr, de - sci, bvec, dipole_mat_k, dipole_mat_kbar)
        
        deallocate (pmat_unocc_occ)

        
        call calc_k_resolved_transition_coefficients(use_ip, dipole_mat_k, dipole_mat_kbar, &
                                           exc_eigvecs, trans_coeffs_k, &
                                           trans_coeffs_kbar)

        ! Sort the transition coffecients by energy if IP approximmation was employed
        if (use_ip) then
            trans_coeffs_k = trans_coeffs_k(ensortidx, :)
            trans_coeffs_kbar = trans_coeffs_kbar(ensortidx, :)
        end if

        n_freqs = input%xs%energywindow%points
        
        allocate (dichroic_tensor_k(3, 3, n_freqs))
        allocate (dichroic_tensor_kbar(3, 3, n_freqs))
        allocate (dichroic_tensor_kkbar(3, 3, n_freqs))

        allocate (osci_strengths_k(3, 3, n_exc))
        allocate (osci_strengths_kbar(3, 3, n_exc))
        allocate (osci_strengths_k_kbar(3, 3, n_exc))

        ! Calculate valley-resolved oscillator strengths
        call calc_k_resolved_oscillator_strengths(input%groundstate%tevecsv, trans_coeffs_k, trans_coeffs_kbar, nk_bse*omega, &
                                                  osci_strengths_k, osci_strengths_kbar, osci_strengths_k_kbar)

        if (mpiglobal%is_root) then
          ! Change sign and divide by pi to make oscilaltor strengths comparable with the ones written in writeoscillator.f90
          call write_k_resolved_osci_strengths(dirnames_oscillator_strengths(1), ham_size, n_exc, nk_bse, -bse_gap, exc_eigvals,-osci_strengths_k / pi)
          call write_k_resolved_osci_strengths(dirnames_oscillator_strengths(2), ham_size, n_exc, nk_bse, -bse_gap, exc_eigvals, -osci_strengths_kbar / pi)
        end if

        ! Generate an evenly spaced frequency grid
        allocate (frequencies(n_freqs))
        call genwgrid(n_freqs, input%xs%energywindow%intv,&
          & input%xs%tddft%acont, 0._dp, w_real=frequencies)

        ! Calculate valley-resolved dichroic tensors
        call calc_k_resolved_dichroic_tensor(frequencies,exc_eigvals,input%xs%broad,osci_strengths_k,osci_strengths_kbar,osci_strengths_k_kbar, dichroic_tensor_k, dichroic_tensor_kbar,dichroic_tensor_kkbar)

        ! write to file
        if (mpiglobal%is_root) then
          call write_dichroic_tensor(input%xs%bse%aresbse,use_ip,input%xs%bse%bsetype,input%xs%screening%screentype,dirnames_dichroic_tensors(1), dichroic_tensor_k, frequencies)
          call write_dichroic_tensor(input%xs%bse%aresbse, use_ip, input%xs%bse%bsetype,input%xs%screening%screentype,dirnames_dichroic_tensors(2), dichroic_tensor_kbar, frequencies)
          call write_dichroic_tensor(input%xs%bse%aresbse,use_ip, input%xs%bse%bsetype,input%xs%screening%screentype,dirnames_dichroic_tensors(3), dichroic_tensor_kkbar, frequencies)

            write (unitout, '(" Dichroic tensors written to file.")')
        end if

    end subroutine

    !> Reads the matrix elements of the momentum operator \( \mathbf{p} \)
    !> between occupied (valence) and unoccupied (conduction) states
    !> for all states considered in BSE.
    !> The momentum matrix elements are defined as (see e.g. Eq. 2)
    !>
    !>\[
    !>      \langle \psi_{v\mathbf{k}}
    !>              | \hat{\mathbf{p}} |
    !>              \psi_{c\mathbf{k}}\rangle
    !>                                          \].
    !>
    !> Matrix elements were written to file in task 'writepmatxs'.
    !> Note: renormalisation works only for materials with a gap.
    subroutine read_momentum_matrix_elements(bse_limits, use_gw, k_vecs_lat, &
                                             eigvals_sv, eigvals_qp, pmat_unocc_occ)

        use m_getpmat, only: getpmat
        use constants, only: zzero, zone, zi

        !> Limits for occupied and unoccpied states considered in BSE
        integer, intent(in) :: bse_limits(:, :)
        !> k-vectors of the BSE calculation in lattice coordinates
        real(dp), intent(in) :: k_vecs_lat(:, :)
        !> Momentum matrix elements between unoccupied and occupied states
        complex(dp), allocatable, intent(out) :: pmat_unocc_occ(:, :, :, :)
        !> True if eigenvalues stem from a GW calculation
        logical, intent(in) :: use_gw
        !> Second variational eigenvalues (needed for renormalisation after GW)
        real(dp), optional, intent(in) :: eigvals_sv(:, :)
        !> GW quasi-particle eigenvalues
        real(dp), optional, intent(in) :: eigvals_qp(:, :)
        !> Running index occupied bands
        integer :: iocc
        !> Running index unoccupied bands
        integer :: iunocc
        !> Absolute index occupied bands
        integer :: ioabs
        !> Absolute index unoccupied bands
        integer ::  iuabs
        !> Running index k-point
        integer :: ik
        integer, parameter :: ik_ref = 1
        !> Number of occupied bands
        integer :: n_occ
        !> Number of unoccupied bands
        integer ::  n_unocc
        !> Absolute index of lowest occupied band
        integer :: iocc_lowest
        !> Absolute index of lowest unoccupied band
        integer :: iunocc_lowest
        !> Absolute index of highest occupied band
        integer ::  iocc_highest
        !> Absolute index of highest unoccupied band
        integer :: iunocc_highest
        !> Number of k-points used in BSE
        integer :: nkpoints_bse

        nkpoints_bse = size(bse_limits, dim=2)

        n_unocc = bse_limits(2, ik_ref) - bse_limits(1, ik_ref) + 1
        n_occ = bse_limits(4, ik_ref) - bse_limits(3, ik_ref) + 1

        allocate (pmat_unocc_occ(3, n_unocc, n_occ, nkpoints_bse))

        do ik = 1, nkpoints_bse

            iunocc_lowest = bse_limits(1, ik)
            iunocc_highest = bse_limits(2, ik)
            iocc_lowest = bse_limits(3, ik)
            iocc_highest = bse_limits(4, ik)

            call getpmat(ik, k_vecs_lat,&
                  & iunocc_lowest, iunocc_highest, &
                  iocc_lowest, iocc_highest,&
                  & .true., 'PMAT_XS.OUT', pmat_unocc_occ(:, :, :, ik))

        end do
        
        if (use_gw) then
         call renormalise_momentum_matrix_elements(pmat_unocc_occ, eigvals_sv, eigvals_qp,iocc_lowest, iocc_highest)
        end if
    end subroutine

    !> Din: Renormalise momentum matrix elements according to Del Sole PRB48, 11789(1993)
    !> \frac{v^        ext{QP}_{okuk}}{E_uk - E_ok} \appprox \frac{p^        ext{LDA}_{okuk}}{e_uk - e_ok}
    !>   In the case that we use the quasi-particle energies E but the LDA eigenfunctions:
    !>   1) evalsv contains the E, while eval0 contains the e.
    !>   2) The BSE diagonal contains \Delta E
    !>   3) But the q->0 limit of <o|exp^{-iqr}|u> is still expressed in terms of
    !>      \frac{p^        ext{LDA}_{okuk}}{e_uk - e_ok}, instead of v and E.
    !>   The following scales the LDA momentum matrix elements such that 3) is
    !>   true for calculations on top of LDA or GW (Eigenvalues only).
    subroutine renormalise_momentum_matrix_elements(pmat, eigvals_sv, eigvals_qp,id_lowest_occ_band,id_lowest_unocc_band)
        
        !> Momentum matrix elements between unoccupied and occupied states
        complex(dp), allocatable, intent(out) :: pmat(:, :, :, :)
        !> Second variational eigenvalues (needed for renormalisation after GW)
        real(dp), optional, intent(in) :: eigvals_sv(:, :)
        !> GW quasi-particle eigenvalues
        real(dp), optional, intent(in) :: eigvals_qp(:, :)
        !> IDs of lowest occupied / unoccupied bands
        integer, intent(in) :: id_lowest_occ_band, id_lowest_unocc_band
        !> Running index occupied bands
        integer :: iocc
        !> Running index unoccupied bands
        integer :: iunocc
        !> Running index k-point
        integer :: ik
        !> Number of occupied bands
        integer :: n_occ
        !> Number of k-points used in BSE
        integer :: nkpt
        !> Number of unoccupied bands
        integer ::  n_unocc

        nkpt = size(pmat,dim=4)
        n_occ = size(pmat,dim=3)
        n_unocc = size(pmat,dim=2)

        do ik = 1, nkpt
        do iocc = 1, n_occ
            do iunocc = 1, n_unocc
                pmat(:, iunocc, iocc, ik) = pmat(:, iunocc, iocc, ik)&
                  &*(eigvals_sv(id_lowest_unocc_band + iunocc - 1, ik) &
                        - eigvals_sv(id_lowest_occ_band + iocc - 1, ik))&
                  &/(eigvals_qp(id_lowest_unocc_band + iunocc - 1, ik) &
                        - eigvals_qp(id_lowest_occ_band + iocc - 1, ik))
            end do
        end do
    end do
    end subroutine 


    !> Generates the resonant Dipole matrix elements \(   \mathbf{D}_{vc\mathbf{k}}  \)
    !> according to Eq. 61 in Vorwerk et al. as
    !>
    !>\[
    !>   \mathbf{D}_{vc\mathbf{k}} = i \frac{\langle \psi_{v\mathbf{k}} | \hat{\mathbf{p}} | \psi_{c\mathbf{k}}\rangle}
    !>                                      {\epsilon_{c\mathbf{k}} - \epsilon_{v\mathbf{k}}}
    !>
    !> Modified from m_setup_dmat.f90
    subroutine calc_k_resolved_dipole_matrix_elements(pmat_unocc_occ, smap, smap_rel, &
                                     k_vecs, energy_diffs, &
                                     reciprocal_lattice, &
                                                      dipole_mat_k, dipole_mat_kbar)
        use constants, only: zzero, zone, zi
        use grid_utils, only: point_in_triangle

        !> K valley-resolved dipole matrix elements
        complex(dp), allocatable, intent(out) :: dipole_mat_k(:, :)
        !> K_bar valley-resolved dipole matrix elements
        complex(dp), allocatable, intent(out) :: dipole_mat_kbar(:, :)
        !> Momentum matrix elements between unoccpied and occupied states
        complex(dp), intent(in) :: pmat_unocc_occ(:, :, :, :)
        !> Map between 1D index of transitions to absolute indices (v_abs,c_abs,k_abs)
        integer, intent(in) :: smap(:, :)
        !> Map between 1D index of transitions to relative indices (v_rel,c_rel,k_rel) (e.g. first used valence band in BSE has index 1)
        integer, intent(in) :: smap_rel(:, :)
        !> k-vectors employed in BSE calulcation in cartesian coordinates
        real(dp), intent(in) :: k_vecs(:, :)
        !> Energy diffs conduction - valence states for all BSE transitions
        real(dp), intent(in)  ::   energy_diffs(:)
        !> Reciprocal lattice vectors
        real(dp), intent(in) :: reciprocal_lattice(3, 3)

        !> Relative index of considered occupied states
        integer :: io
        !> Relative index of considered unoccupied states
        integer ::  iu
        !> Absolute index of considered unoccupied states
        integer :: iuabs
        !> Absolute index of considered occupied states
        integer :: ioabs
        !> Index of reduced k-point
        integer :: ik
        !> Index of non-reduced k-point
        integer :: iknr
        !> Running index transitions
        integer :: a1
        !> Cartesian coordinates of k_vector
        real(dp) :: vkc0(3)
        !> Size of BSE Hamiltonian (ham_size x ham_size)
        integer :: ham_size
        !> One dipole matrix element, will be assigned to the valleys
        complex(dp) :: dipole_element(3)

        ham_size = size(smap, dim=2)

        allocate (dipole_mat_k(ham_size, 3))
        allocate (dipole_mat_kbar(ham_size, 3))

        !$omp PARALLEL DO &
        !$omp& DEFAULT(SHARED), PRIVATE(a1, iknr, iu, io, ik, vkc0, dipole_element)
        do a1 = 1, ham_size

            ! Absolute indices
            iknr = smap(3, a1)

            ! Relative indices
            iu = smap_rel(1, a1)
            io = smap_rel(2, a1)
            ik = smap_rel(3, a1)

            ! Note: The scissor does not enter here, so subtract it again.
            dipole_element = zi*pmat_unocc_occ(:, iu, io, ik)/ energy_diffs(a1) 

            ! Get cartesian coordinates of current k-point
            ! and check if it's in valley around high-symmetry point K
            if (point_in_triangle(reciprocal_lattice(1:2, 1), &
                                    reciprocal_lattice(1:2, 2), &
                                    k_vecs(1:2, iknr))) then
                dipole_mat_k(a1, :) = dipole_element

                dipole_mat_kbar(a1, :) = zzero
            else
                dipole_mat_kbar(a1, :) = dipole_element

                dipole_mat_k(a1, :) = zzero

            end if

        end do
        !$omp END PARALLEL DO

    end subroutine calc_k_resolved_dipole_matrix_elements

    !> Calculates the valley-resolved transition coefficients \( \mathbf{t}^\lambda \)
    !> (see Eq. (2) of Caruso et al. and Eq. (62) of Vorwerk et al.)
    !>
    !> \[ \mathbf{t}^\lambda = -i \sum_{vc}\sum^{\mathcal{P}_{\rm K}}_\mathbf{k} [A^\lambda_{vc\mathbf{k}}] ^*\mathbf{D}_{vc\mathbf{k}}^*
    !>
    !>          \]
    !>
    !> with the BSE eigenvectors \(  A^\lambda_{vc\mathbf{k}} \) and the dipole elements
    !> \(   \mathbf{D}_{vc\mathbf{k}}  \) (see routine calc_k_resolved_dipole_matrix_elements for definition.)
    subroutine calc_k_resolved_transition_coefficients(use_ip, dipole_mat_k, &
                                                       dipole_mat_kbar, exc_eigvecs, &
                                             trans_coeff_k, trans_coeff_kbar)

        use xlapack, only: matrix_multiply

        implicit none

        !> BSE eigenvectors
        complex(dp), intent(in) :: exc_eigvecs(:, :)
        !> K valley-resolved oscillator strength
        complex(dp), allocatable, intent(out) :: trans_coeff_k(:, :)
        !> K_bar valley-resolved oscillator strength
        complex(dp), allocatable, intent(out) :: trans_coeff_kbar(:, :)
        !> K valley-resolved dipole matrix element
        complex(dp), intent(in) :: dipole_mat_k(:, :)
        !> K_bar valley-resolved dipole matrix element
        complex(dp), intent(in) :: dipole_mat_kbar(:, :)

        !> True if independent-particle approximation was employed in BSE
        logical, intent(in) :: use_ip
        !> Number of solution of BSE diagonalisation
        integer ::  n_exc

        n_exc = size(exc_eigvecs, dim=1)

        allocate (trans_coeff_kbar(n_exc, 3))
        allocate (trans_coeff_k(n_exc, 3))

        ! BSE eigenvectors are diagonal ( i.e. \delta^\lambda_{vc\mathbf{k}} ) in IP approximation
        if (use_ip) then

            trans_coeff_k =  conjg(zi*dipole_mat_k)
            trans_coeff_kbar =  conjg(zi*dipole_mat_kbar)

        else
            call matrix_multiply(exc_eigvecs,  conjg(zi*dipole_mat_k), trans_coeff_k, trans_A='C')
            call matrix_multiply(exc_eigvecs,   conjg(zi*dipole_mat_kbar), trans_coeff_kbar, trans_A='C')
        end if

    end subroutine calc_k_resolved_transition_coefficients

    !> Calculates the valley-resolved excitonic oscillator strengths defined as
    !>(see Eq. 3 and following text of Caruso et al.)
    !>
    !> \[
    !>        \frac{4\pi^2}{V} (t^{\lambda,\alpha}_{\rm K}) ^*t^{\lambda, \beta}_{\rm K}
    !>                  \]
    !>
    !> with the valley-resolved transition coefficients \( t^{\lambda,\alpha}_{\rm K}\) .
    !> By making the replacements (K -> K_bar, K -> K_bar) and (K -> K, K -> K_bar)
    !> we can obtain the oscillator strengths for the K_bar valley and the "mixed" oscillator strength.
    subroutine calc_k_resolved_oscillator_strengths(use_spinorbit,  trans_coeffs_k, trans_coeffs_kbar, crystal_volume, osci_strengths_k, osci_strengths_kbar, osci_strengths_k_kbar)
        
        use constants, only: zone, zi, pi
        use modmpi

        implicit none

        !> True if spin-orbit coupling was employed for ground state
        logical, intent(in) :: use_spinorbit

        !> Crystal volume (V = N_kpt * Omega_{unit-cell})
        real(dp), intent(in) :: crystal_volume
        !> K-bar-valley resolved transition coeffient
        complex(dp), intent(in) :: trans_coeffs_k(:, :)
        !> K-bar-valley resolved transition coeffient
        complex(dp), intent(in) :: trans_coeffs_kbar(:, :)
        !> K-valley resolved oscillator strength
        complex(dp), intent(out) :: osci_strengths_k(:, :, :)
        !> K-valley resolved oscillator strength
        complex(dp), intent(out) ::     osci_strengths_kbar(:, :, :)
        !> Mixed valley resolved oscillator strength
        complex(dp), intent(out) ::     osci_strengths_k_kbar(:, :, :)

        !> Number of excitons found from diagonalizing the BSE Hamiltonian
        integer :: n_exc
        !> Momentum transfer index
        integer, parameter :: iqmt = 1
        !> Running index excitons
        integer :: iexc
        !> Running indices optical cartesian directions
        integer :: jopt, iopt
        !> Running index frequencies
        integer ::ifreq
        !> Prefactor pref = 4*pi / V (times 2 for no spin splitting) 
        real(dp) :: pref  
        
        n_exc = size(trans_coeffs_k, dim=1)

        ! Adjusting prefactor, the factor 2 accounts for spin degeneracy
        pref = -4._dp*pi**2/crystal_volume
        if (.not. use_spinorbit)  pref = pref*2._dp
           

        do iopt = 1, 3
            do jopt = 1, 3

                osci_strengths_k(iopt, jopt, :) = pref*conjg(trans_coeffs_k(:, iopt))*trans_coeffs_k(:, jopt)
                osci_strengths_kbar(iopt, jopt, :) = pref*conjg(trans_coeffs_kbar(:, iopt))*trans_coeffs_kbar(:, jopt)
                osci_strengths_k_kbar(iopt, jopt, :) = pref*conjg(trans_coeffs_k(:, iopt))*trans_coeffs_kbar(:, jopt)

            end do
        end do

    end subroutine calc_k_resolved_oscillator_strengths

    !> Calculates the valley-resolved dichroic tensor \(  \Xi^{\rm K}_{\alpha\beta}    \) according to Eq. (3) of Caruso et al.:
    !>
    !> \[
    !>
    !>      \xi^{\rm K}_{\alpha\beta}  = \frac{4\pi^2}{V} 
    !>                                   \sum_\lambda (t^{\lambda,\alpha}_{\rm K})^* t^{\lambda, \beta}_{\rm K} \delta(E^\lambda - \omega)
    !>
    !>    \]
    !> 
    !> with the excitonic oscillator strengths \( \frac{4\pi^2}{V} (t^{\lambda,\alpha}_{\rm K}) ^*t^{\lambda, \beta}_{\rm K} \) and
    !> eigenenergies \( E^\lambda \). Note that in exciting the frequency term \( \delta(E^\lambda - \omega) \) is replaced by a Lorentzian
    !> with poles at both \( \pm E^\lambda \) in order to account for resonant and anti-resonant part of the BSE spectrum (see Vorwerk et al. for a discussion on this.)
    subroutine calc_k_resolved_dichroic_tensor(frequencies, exc_eigvals, broadening,osci_strengths_k,osci_strengths_kbar, osci_strengths_k_kbar, dichroic_tensor_k, dichroic_tensor_kbar, dichroic_tensor_kkbar)
        use tensor_contractions, only: complex_tensor_contraction_dp
        use distributions, only: lorentzian
        
        !> K-bar valley resolved dichroic tensor
        complex(dp), intent(in) ::  osci_strengths_k(:, :, :)
        !> K-bar valley resolved dichroic tensor
        complex(dp), intent(in) :: osci_strengths_kbar(:, :, :)
        !> Mixed valley resolved dichroic tensor
        complex(dp), intent(in) ::   osci_strengths_k_kbar(:, :, :)
        
        !> frequencies for generation of dichroic tensors
        real(dp), intent(in) :: frequencies(:)
        !> Exciton eigenenergies
        real(dp), intent(in) :: exc_eigvals(:)
   
        !> Broadening parameter for delta function
        real(dp), intent(in) :: broadening
        !> K-bar valley resolved dichroic tensor
        complex(dp), intent(out) :: dichroic_tensor_kbar(:, :, :)
        !> K-K-bar mixed dichroic tensor
        complex(dp), intent(out) ::  dichroic_tensor_kkbar(:, :, :)
                !> K-valley resolved dichroic tensor
        complex(dp), intent(out) :: dichroic_tensor_k(:, :, :)
        complex(dp), allocatable :: energy_denominator(:, :)
         !> Number of solution of BSE diagonalisation
        integer ::  n_exc
        !> Number of frequencies in frequency grid
        integer ::  n_freqs
        !> Running index excitons
        integer :: iexc


        ! enw_{w, \lambda} = 1/(w - E_\lambda + i\delta) + 1/(-w - E_\lambda - i\delta)
        n_exc = size(osci_strengths_k, dim=3)
        n_freqs = size(frequencies)
        allocate (energy_denominator(n_exc, n_freqs))

        do iexc = 1, n_exc
            energy_denominator(iexc, :) = &
                -lorentzian(broadening, frequencies, exc_eigvals(iexc)) &
                + lorentzian(broadening, frequencies, -exc_eigvals(iexc))
        end do

        call complex_tensor_contraction_dp(osci_strengths_k, &
                                           shape(osci_strengths_k), &
                                           energy_denominator, shape(energy_denominator), &
                                           dichroic_tensor_k, shape(dichroic_tensor_k))
        call complex_tensor_contraction_dp(osci_strengths_kbar, &
                                           shape(osci_strengths_kbar), &
                                           energy_denominator, shape(energy_denominator), &
                                           dichroic_tensor_kbar, shape(dichroic_tensor_kbar))
        call complex_tensor_contraction_dp(osci_strengths_k_kbar, &
                                           shape(osci_strengths_k_kbar), &
                                           energy_denominator, shape(energy_denominator), &
                                           dichroic_tensor_kkbar, shape(dichroic_tensor_kkbar))

    end subroutine



    !> Writes the dichroic tensor to human-readable file.
    subroutine write_dichroic_tensor(use_antiresonant, use_ip, bse_type, screen_type, dirname, dichroic_tensor, frequencies)
        use modmpi
        use m_genfilname

        use m_writeeps

         !> True if independent-particle approximation was employed in BSE
        logical, intent(in) :: use_ip
        !> True if anti-resonant part was used in BSE
        logical, intent(in) :: use_antiresonant
        !> Frequency grid of the spectrum
        real(dp), intent(in) :: frequencies(:)
        !> Dichroic tensor
        complex(dp), intent(in) :: dichroic_tensor(:, :, :)
        !> Directory name where to store the file
        character(256), intent(in) :: dirname
        !> BSE type (Tamm-Dancoff or not etc.)
        character(256), intent(in) :: bse_type
        !> Screening type (full etc.)
        character(256), intent(in) :: screen_type
        !> Filename
        character(256) :: filename
        !> Strings depending on BSE and screening details
        character(256) :: tdastring, bsetypestring, scrtypestring
        !> Number of frequencies in frequency grid
        integer :: n_freqs
        !> Running indices optical indices (=cartesian directions)
        integer :: jopt, iopt
        
        ! Currently VSCD is only implemented for q=0 independent-particle
        ! approximation or in the Tamm-Dancoff approximation (TDA)
        if (use_ip) then
            tdastring = '' 
        else 
            tdastring = "-TDA-BAR"
        end if

        bsetypestring = '-'//trim(bse_type)//trim(tdastring)
        scrtypestring = '-'//trim(screen_type)

        n_freqs = size(frequencies)

        do iopt = 1, 3
            do jopt = 1, 3

              ! Generate File names for resulting quantities
              call genfilname(basename=trim(dirname), tq0=.true., oc1=iopt, oc2=jopt,&
                & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
                    & nar=.not. use_antiresonant, filnam=filename)

              ! Write optical functions to file
                call writeeps(1, iopt, jopt, frequencies, dichroic_tensor(iopt, jopt, :), trim(dirname)//'/'//trim(filename))
            end do
        end do

    end subroutine write_dichroic_tensor

    !> Writes the valley-resolved oscillator strengths to human-readable file.
    subroutine write_k_resolved_osci_strengths(dirname, hamsize, nexc, nk, band_gap, exciton_eigvals, osci_strengths)
       use modmpi
        use modinput, only: input
        use modxs, only: bsed, escale, ivgmt, vqlmt, vgcmt, vqcmt
        use m_genfilname
        use m_getunit
        use m_write_hdf5
        use modxs, only: symt2, ivgigq, sptclg
        use constants, only: zzero
        use mod_lattice, only: omega

        !> BSE Hamiltonian size
        integer, intent(in) :: hamsize
        !> Number of excitations found in BSE (mostly equal to hamiltonian size)
        integer, intent(in) :: nexc
        !> Number of k-points used in BSE 
        integer, intent(in) :: nk
        !> Band gap (including scissor correction if applied)
        real(dp), intent(in) :: band_gap
        !> BSE eigenvalues
        real(dp), intent(in) :: exciton_eigvals(hamsize)
        !> Name of directory where to store files 
        character(256), intent(in) :: dirname

        !> Running indices optical coordinates (= cartesian directions)
        integer :: iopt, jopt
        !> Running index excitons
        integer ::  iexc
        !> File unit
        integer :: unexc
        integer, parameter ::  iq_gamma = 1
      character(256) :: fnexc, frmt, tdastring, bsetypestring, scrtypestring
      character(256) :: syscommand

        complex(dp), intent(in) :: osci_strengths(3, 3, nexc)
        complex(dp) :: osci_strengths_symm(3, 3, nexc)

      ! Generate file name
        if (input%xs%bse%coupling) then
            tdastring = ''
      else
            if (input%xs%bse%chibarq) then
                tdastring = "-TDA-BAR"
        else
                tdastring = "-TDA"
        end if
      end if
        if (input%xs%bse%bsetype == "IP") then
            tdastring = ''
      end if
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)


        
        ! symmetrize the oscillator strength
        do iopt = 1, 3
            do jopt = 1, 3
            if (.NOT. input%xs%BSE%chibar0) then 
                osci_strengths_symm(iopt, jopt, :) = osci_strengths(iopt, jopt, :)
            else
                call symt2app(iopt, jopt, nexc, symt2, osci_strengths, osci_strengths_symm(iopt, jopt, :))
            end if

          
        call genfilname(basename=dirname, tq0=.true., oc1=iopt, oc2=jopt,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
                        & nar=.not. input%xs%bse%aresbse, filnam=fnexc)
        
            fnexc = trim(dirname)//'/'//trim(fnexc)

        ! Write out exciton energies and oscillator strengths
        call getunit(unexc)
            open (unexc, file=fnexc, form='formatted', action='write', status='replace')
            write (unexc, '("#",1x,"Excitonic eigen energies and oscillator strengths")')
            write (unexc, '("#")')
            write (unexc, '("# Momentum transfer Q=G+q in lattice cooridnates")')
            write (unexc, '("# G:",3i4)') ivgmt(1:3, iq_gamma)
            write (unexc, '("# q:",3f12.7)') vqlmt(1:3, iq_gamma)
            write (unexc, '("# Momentum transfer Q=G+q in Cartesian cooridnates")')
            write (unexc, '("# G:",3f12.7)') vgcmt(1:3, iq_gamma)
            write (unexc, '("# q:",3f12.7)') vqcmt(1:3, iq_gamma)
            write (unexc, '("# Norm2(G+q)",f12.7)') norm2(vgcmt(:, iq_gamma) + vqcmt(:, iq_gamma))
            write (unexc, '("#")')
            write (unexc, '("# Energy scale",f12.7)') escale
            write (unexc, '("# E_shift : ", SP, E23.16)') band_gap*escale
            write (unexc, '("#")')

            frmt = '(a1,a7,5(1x,a23))'
            write (unexc, frmt) "#", "Nr.",&
          & "E",&
          & "E+E_shift",&
          & "|Osc. Str.|",&
          & "Re(Osc. Str. Res.)",&
          & "Im(Osc. Str. Res.)"
        !frmt='(I8,5(1x,E23.16))'
            frmt = '(I8,2(1x,F23.16),3(1x,E23.16))'
            do iexc = 1, nexc
                write (unexc, frmt) iexc,&
              & exciton_eigvals(iexc)*escale,&
                        & (exciton_eigvals(iexc) + band_gap)*escale,&
                        & abs(osci_strengths_symm(iopt, jopt, iexc)),&
                        & dble(osci_strengths_symm(iopt, jopt, iexc)),&
                        & aimag(osci_strengths_symm(iopt, jopt, iexc))
        end do

            close (unexc)

            end do 
        end do 

    end subroutine write_k_resolved_osci_strengths

end module dichroism
