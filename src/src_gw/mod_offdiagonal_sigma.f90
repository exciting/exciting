!> Module providing all equations to allocate
!> save, delete, and compute the offdiagonal
!> elements.
module mod_offdiagonal_selfenergy

    use precision,  only: i32, dp, str_512, long_int
    use constants,  only: pi, zzero, real_zero
    use iso_c_binding, only: c_loc 
#include "offload.fpp"
#ifdef USEOMP
    use omp_lib
#endif
    
    implicit none
    private

    !----------------------------------------------------------------------
    !> Root filename for the exchange off-diagonal self-energy output
    character(len=*), private, parameter :: exchange_offdiagonal_selfenergy_rootname = 'SIGMAX_OFFDIAGONAL_K'

    !> Root filename for the correlation off-diagonal self-energy output
    character(len=*), private, parameter :: correlation_offdiagonal_selfenergy_rootname = 'SIGMAC_OFFDIAGONAL_K'

    !----------------------------------------------------------------------
    !> Off-diagonal terms of the exchange (static) self-energy
    complex(dp), private, allocatable, target :: sigma_exchange_offdiagonal(:,:)

    !> Off-diagonal terms of the correlation (dynamic) self-energy
    complex(dp), private, allocatable, target :: sigma_correlation_offdiagonal(:,:,:)

    ! Common functions
    public :: delete_offdiagonal_selfenergy
    public :: add_offdiag_selfenergy_at_energies_to_optimized_xc_potential

    ! Functions related to the correlation contribution
    public :: init_offdiagonal_selfenergy_correlation
    public :: mpi_reduce_offdiagonal_selfenergy_correlation
    public :: write_offdiagonal_selfenergy_correlation
    public :: read_offdiagonal_selfenergy_correlation
    public :: add_q_omega_contrib_to_offdiagonal_selfenergy_correl_at_ik
    

    ! Functions related to the exchange contribution
    public :: init_offdiagonal_selfenergy_exchange
    public :: mpi_reduce_offdiagonal_selfenergy_exchange
    public :: write_offdiagonal_selfenergy_exchange
    public :: read_offdiagonal_selfenergy_exchange
    public :: add_q_contrib_to_offdiagonal_selfenergy_exchange_at_ik

contains

    !> Deallocate all off-diagonal self-energy arrays.
    !> Safe cleanup of allocated tensors.
    subroutine delete_offdiagonal_selfenergy()
        if (allocated(sigma_exchange_offdiagonal))     deallocate(sigma_exchange_offdiagonal)
        if (allocated(sigma_correlation_offdiagonal)) then
            OMP_OFFLOAD target exit data map(always, delete: sigma_correlation_offdiagonal)
            deallocate(sigma_correlation_offdiagonal)
        end if
    end subroutine delete_offdiagonal_selfenergy


    !> Initialize exchange (static) off-diagonal self-energy.
    !> Allocates and zeroes the exchange self-energy tensor.
    subroutine init_offdiagonal_selfenergy_exchange(init_kpt, final_kpt)
        use mod_eigenvalue_occupancy, only: nstfv
        use math_utils, only: number_of_upper_triangle_elements

        !> Starting k-point index
        integer(i32), intent(in) :: init_kpt
        !> Final k-point index
        integer(i32), intent(in) :: final_kpt

        ! Deallocate any previously allocated exchange array
        if (allocated(sigma_exchange_offdiagonal)) deallocate(sigma_exchange_offdiagonal)

        ! Allocate and initialize to zero
        allocate(sigma_exchange_offdiagonal(number_of_upper_triangle_elements(nstfv), init_kpt:final_kpt), source=zzero)
    end subroutine init_offdiagonal_selfenergy_exchange


    !> Initialize correlation (dynamic) off-diagonal self-energy.
    !> Allocates and zeroes the correlation self-energy tensors.
    subroutine init_offdiagonal_selfenergy_correlation(init_kpt, final_kpt, init_freq, final_freq)
        use mod_eigenvalue_occupancy, only: nstfv
        use math_utils, only: number_of_upper_triangle_elements
        use mod_device_offload, only: device_world

        !> Starting k-point index
        integer(i32), intent(in) :: init_kpt
        !> Final k-point index
        integer(i32), intent(in) :: final_kpt
        !> Starting frequency index
        integer(i32), intent(in) :: init_freq
        !> Final frequency index
        integer(i32), intent(in) :: final_freq

        integer(i32) :: device_id 

        device_id = device_world%get_device()

        ! Deallocate any previously allocated correlation arrays
        if (allocated(sigma_correlation_offdiagonal)) then
            OMP_OFFLOAD target exit data map(delete: sigma_correlation_offdiagonal) if(omp_target_is_present(c_loc(sigma_correlation_offdiagonal), device_id) /= 0)
            deallocate(sigma_correlation_offdiagonal)
        end if

        ! Allocate arrays and initialize to zero
        allocate(sigma_correlation_offdiagonal(number_of_upper_triangle_elements(nstfv), init_freq:final_freq, init_kpt:final_kpt), source=zzero)
        OMP_OFFLOAD target enter data map(always, to: sigma_correlation_offdiagonal)

    end subroutine init_offdiagonal_selfenergy_correlation
    
    !> MPI reduces sigma_correlation_offdiagonal for the given MPI communicator
    !> Note that after the reduction for device aware calculations, the device
    !> memory is not syncronized.
    subroutine mpi_reduce_offdiagonal_selfenergy_correlation(mpi_environment)
        use mod_mpi_gw,   only: mpi_sum_array
        use exciting_mpi, only: mpiinfo
        !> MPI environment
        type(mpiinfo), intent(in) :: mpi_environment

        OMP_OFFLOAD target update from(sigma_correlation_offdiagonal)
        call mpi_sum_array(sigma_correlation_offdiagonal, mpi_environment, .false.)

    end subroutine mpi_reduce_offdiagonal_selfenergy_correlation

    !> Dumps the offdiagonal selfenergy correlation in a file in the given format
    subroutine write_offdiagonal_selfenergy_correlation( ik, file_format )
        use gw_io, only: build_file_name, write_to_file
        !> Index of the current k-point
        integer(i32), intent(in) :: ik
        !> Format of the file where to print. It can be e.g. 'text' or 'binary'
        character(len=*), intent(in) :: file_format

        character(len=str_512) :: file_name
        integer(i32) :: lbounds(3)

        ! Use the correlation rootname instead of exchange
        call build_file_name( correlation_offdiagonal_selfenergy_rootname, ik, file_name )
       
        lbounds = lbound(sigma_correlation_offdiagonal)

        ! Pass the correlation array instead of the exchange array
        call write_to_file( file_name, sigma_correlation_offdiagonal(:, :, ik), &
                            lbounds(1:2), file_format)

    end subroutine write_offdiagonal_selfenergy_correlation

    !> MPI reduces sigma_exchange_offdiagonal for the given MPI communicator
    subroutine mpi_reduce_offdiagonal_selfenergy_exchange(mpi_environment)
        use mod_mpi_gw,   only: mpi_sum_array
        use exciting_mpi, only: mpiinfo
        !> MPI environment
        type(mpiinfo), intent(in) :: mpi_environment

        call mpi_sum_array( sigma_exchange_offdiagonal, mpi_environment, .false.)

    end subroutine mpi_reduce_offdiagonal_selfenergy_exchange

    !> Dumps the offdiagonal selfenergy exchange in a file in the given format
    subroutine write_offdiagonal_selfenergy_exchange( ik,  file_format)
        use gw_io, only: build_file_name, write_to_file
        !> Index of the current k-point
        integer(i32), intent(in) :: ik
        !> Format of the file where to print. It can be e.g. 'text' or 'binary'
        character(len=*), intent(in) :: file_format

        character(len=str_512) :: file_name

        call build_file_name( exchange_offdiagonal_selfenergy_rootname, ik, file_name )
        call write_to_file( file_name, sigma_exchange_offdiagonal(:, ik), lbound( sigma_exchange_offdiagonal, dim=1), file_format)

    end subroutine write_offdiagonal_selfenergy_exchange

    !> Computes the q,k contribution to the off-diagonal exchange self-energy at k-point ik
    subroutine add_q_contrib_to_offdiagonal_selfenergy_exchange_at_ik( &
                                            ikp, ispace_init, ispace_final, mdim, nomax, jk, &
                                            degenerate_subspaces, minmmat, kiw, corind, idxas, ciw)

        use math_utils,               only: flatten_idx_for_elements_upper_triangle
        use mod_eigenvalue_occupancy, only: nstfv
    
        implicit none
        
        !> Index of the k-point where the exchange self-energy is evaluated
        integer(i32), intent(in) :: ikp
        !> Range of degenerate subspaces to be considered at k-point ikp
        integer(i32), intent(in) :: ispace_init, ispace_final
        !> Total number of states (valence + core) entering the exchange sum
        integer(i32), intent(in) :: mdim
        !> Number of occupied valence states; states > nomax correspond to core levels
        integer(i32), intent(in) :: nomax
        !> Index of the k-point used in the Brillouin-zone summation
        integer(i32), intent(in) :: jk
        !> Degeneracy information for bands at each k-point:
        !> (1, ispace, ikp) → lowest band index of the subspace
        !> (2, ispace, ikp) → highest band index of the subspace
        !> (3, ispace, ikp) → size (degeneracy) of the subspace
        integer(i32), intent(in) :: degenerate_subspaces(:, :, :)
        !> Expansion coefficients at k,q
        !> Indexed as (basis, band_i, band_j)
        complex(dp), allocatable, intent(in) :: minmmat(:, :, :)
        !> k-point weights / occupation factors for valence states
        !> kiw(ie, jk) gives the weight of valence state ie at k-point jk
        real(dp), intent(in) :: kiw(:, :)
        !> Mapping between core-state index and atomic quantum numbers
        !> corind(icg,1) → atomic species index
        !> corind(icg,2) → atom index within species
        !> corind(icg,6) → core level index
        integer(i32), intent(in) :: corind(:, :)
        !> Mapping from (atom, species) to a global atom-species index
        integer(i32), intent(in) :: idxas(:, :)
        !> Occupation weights for core states
        !> ciw(ic, ias) gives the weight of core level ic on atom-species ias
        real(dp), intent(in) :: ciw(:, :)

        integer(i32) :: ispace, jspace
        integer(i32) :: lowband, upband, size_deg
        integer(i32) :: lowbandp, upbandp, size_degp
        integer(i32) :: first_band_idx
        integer(i32) :: ie1, ie2, ie3, flatten_idx
        integer(i32) :: icg, is, ia, ic, ias
        complex(dp) :: sx, mvm

        first_band_idx = degenerate_subspaces(1, 1, ikp)
        
        !$omp parallel do default(none) schedule(dynamic) &
        !$omp private(ispace, jspace, lowband, upband, size_deg, lowbandp, upbandp, size_degp) &
        !$omp private(ie1, ie2, ie3, sx, mvm, icg, is, ia, ic, ias, flatten_idx) &
        !$omp shared(sigma_exchange_offdiagonal, degenerate_subspaces, minmmat) &
        !$omp shared(kiw, corind, idxas, ciw, ikp, nomax, mdim, jk) &
        !$omp shared(ispace_init, ispace_final, first_band_idx, nstfv)
        do ispace = ispace_init, ispace_final
            do jspace = ispace_init, ispace_final
        
            lowband   = degenerate_subspaces(1, ispace, ikp)
            upband    = degenerate_subspaces(2, ispace, ikp)
            size_deg  = degenerate_subspaces(3, ispace, ikp)
        
            lowbandp  = degenerate_subspaces(1, jspace, ikp)
            upbandp   = degenerate_subspaces(2, jspace, ikp)
            size_degp = degenerate_subspaces(3, jspace, ikp)
        
            ! At a given k-point, degenerate bands span representations of the little group.
            ! Each set of bands with a given degeneracy corresponds to an irreducible
            ! representation (irrep) of the little group.
            !
            ! Subspaces with different degeneracies must belong to different irreps.
            ! By Schur’s lemma, any operator that commutes with all symmetries of the k-point
            ! (such as the Hamiltonian or self-energy Σ) has vanishing matrix elements
            ! between states belonging to different irreps.
            !
            ! Therefore, for states |i> and |j> belonging to subspaces of different
            ! degeneracies (and thus different irreps):
            !     <i | Σ | j> = 0
            !
            ! If two subspaces have the same degeneracy, they may belong to the same irrep
            ! or to distinct irreps of equal dimension; in this case, the matrix elements
            ! between them are not guaranteed to vanish.
            !
            ! While, one should explicitly identify and use the irreducible
            ! representations at each k-point. This is not an easy task
            ! therefore, we use the degeneracy size as a first filter.
            if (size_deg /= size_degp) cycle

            do ie2 = lowbandp, upbandp
                do ie1 = lowband, upband

                ! The flatten index in the upper diagonal, if not there cycle
                flatten_idx = flatten_idx_for_elements_upper_triangle(nstfv, ie1, ie2, first_band_idx, first_band_idx)
                if (flatten_idx < 0) cycle

                sx = zzero
        
                do ie3 = 1, mdim
                    if (ie3 <= nomax) then
                        ! Valence
                        mvm = dot_product(minmmat(:,ie1,ie3), &
                                            minmmat(:,ie2,ie3))
                        sx = sx - kiw(ie3,jk) * mvm
                    else
                        ! Core
                        icg = ie3 - nomax
                        is  = corind(icg,1)
                        ia  = corind(icg,2)
                        ic  = corind(icg,3)
                        ias = idxas(ia,is)
            
                        mvm = dot_product(minmmat(:,ie1,ie3), &
                                            minmmat(:,ie2,ie3))
                        sx = sx - ciw(ic,ias) * mvm
                    end if
                end do

                sigma_exchange_offdiagonal(flatten_idx,ikp) = &
                    sigma_exchange_offdiagonal(flatten_idx,ikp) + sx
        
                end do
            end do
        
            end do
        end do
        !$omp end parallel do
    
    end subroutine add_q_contrib_to_offdiagonal_selfenergy_exchange_at_ik
    
    !> Computes the contribution of the M(k,q)W(q)M(k,q) product to the offdiagonal
    !> correlation self-energy. 
    !>
    !> Note that to reduce memory usage the convolution is done here at this step for each element
    !> instead of over all elements and frequencies, as it would be too much
    !>
    !> \Sigma_{nl\mathbf{k}}^{c}(i\omega) = \frac{1}{N_{c}} \frac{1}{2\pi} \sum_{\mathbf{q}}^{BZ}
    !> \sum_{m} \sum_{i,j} [M_{nm}^{i}(\mathbf{k}, \mathbf{q})]^* ! \int_{0}^{\infty} \frac{2(\epsilon_{m,\mathbf{k}-\mathbf{q}} - i\omega)
    !> [W_{i,j}^{c}(\mathbf{q}, i\omega') - W_{i,j}^{c}(\mathbf{q}, i\omega)]}
    !> {(i\omega - \epsilon_{m,\mathbf{k}-\mathbf{q}})^2 + \omega'^2} d\omega' M_{lm}^{j}(\mathbf{k}, \mathbf{q})
    !> + \frac{1}{2} sgn(\epsilon_{m,\mathbf{k}-\mathbf{q}}) W_{i,j}^{c}(\mathbf{q}, i\omega)
    subroutine add_q_omega_contrib_to_offdiagonal_selfenergy_correl_at_ik( &
                nstart, nend, mstart, mend, nstse, &
                evalfv, evalcr, corind, idxas, efermi, &
                minm, wm, wkq, freq, freq_selfc, &
                nomeg, iom, ikp, jk)

        use mod_frequency,       only: frequency
        use mod_gw_degeneracies, only: band_degeneracy
        use math_utils,          only: flattened_upper_triangle_index_to_element_indexes, &
                                       number_of_upper_triangle_elements
        use modmpi,              only: terminate_if_false
        use modgw,               only: mbsiz, kset

        implicit none

        !> Starting index of external bands (n/l)
        integer(i32),  intent(in) :: nstart
        !> Ending index of external bands (n/l)
        integer(i32),  intent(in) :: nend
        !> Starting index of intermediate bands (m)
        integer(i32),  intent(in) :: mstart
        !> Ending index of intermediate bands (m)
        integer(i32),  intent(in) :: mend
        !> Number of bands in the calculation
        integer(i32),  intent(in) :: nstse
        !> Frequency index count for the convolution
        integer(i32),  intent(in) :: nomeg
        !> Index of the external frequency (iω)
        integer(i32),  intent(in) :: iom
        !> Index of the k-point (irreducible)
        integer(i32),  intent(in) :: ikp
        !> Index of the k-point used in the Brillouin-zone summation (full BZ)
        integer(i32), intent(in)  :: jk

        !> Fermi energy
        real(dp), intent(in) :: efermi

        !> Valence/conduction band energies: evalfv(m, k)
        real(dp), intent(in) :: evalfv(:,:)
        !> Core state energies: evalcr(core_index, state_index)
        real(dp), intent(in) :: evalcr(:,:)
        !> Mapping of core indices: corind(core_group, :) = (is, ia, ic)
        integer(i32),  intent(in) :: corind(:,:)
        !> Atomic state index mapping: idxas(atom indx within species, species idx)
        integer(i32),  intent(in) :: idxas(:,:)

        !> Matrix elements [M^i_{nm}] for the off-diagonal computation
        complex(dp), intent(in) :: minm(mbsiz,nstart:nend,mstart:mend)
        !> Matrix elements [W(q) M^j_{lm}] for the off-diagonal computation. 
        !> If compiled with device offload support, it is a host-residing pointer 
        !> to a device memory
        complex(dp), pointer, intent(in) :: wm(:,:,:)
        !> Pre-factor including the BZ integration weight
        real(dp), intent(in) :: wkq

        !> Frequency data structure (includes freqs and weights womeg)
        type(frequency), intent(in) :: freq
        !> Self-energy frequencies data structure (freqs)
        type(frequency), intent(in) :: freq_selfc


        ! Local variables
        integer(i32)      :: ie1, ie2, ie3, iom2, nstates, jkp
        integer(long_int) :: ntri, itri
        integer(i32)      :: icg, is, ia, ic, ias
        real(dp)          :: energy_ie3 
        complex(dp)       :: mwm_offdiagonal, convolution_contribution, convolution_contribution_p
        complex(dp)       :: energy_factor, z_energy, energy_factor_p, z_energy_p

        nstates =  nend - nstart + 1
        ntri    =  number_of_upper_triangle_elements(nstates)
        jkp     =  kset%ik2ikp(jk)

        OMP_OFFLOAD target enter data map(always, to: freq, freq_selfc)
#if defined(FLANG_OPENMP_DERIVED_TYPE_MAP_BUG_WORKAROUND)
        OMP_OFFLOAD target enter data map(always, to: freq%freqs, freq_selfc%freqs)
#endif  

        ! Loop over the external indices 'n' and 'l' (left and right indices of Sigma)
        OMP_OFFLOAD target has_device_addr(wm)
        !$omp teams distribute parallel do default(none) &
        !$omp shared(nstart, nend, mstart, mend, nstse) &
        !$omp shared(evalfv, evalcr, corind, idxas, efermi) &
        !$omp shared(minm, wm, wkq, freq, freq_selfc, band_degeneracy) &
        !$omp shared(sigma_correlation_offdiagonal, iom, ikp, nomeg, jkp, ntri, nstates) &
        !$omp private(ie1, ie2, ie3, iom2, energy_ie3, z_energy) &
        !$omp private(mwm_offdiagonal, energy_factor, convolution_contribution) &
        !$omp private(icg, is, ia, ic, ias, itri, energy_factor_p, z_energy_p) &
        !$omp private(convolution_contribution_p)
        do itri = 1, ntri

            ! Compute the global indexes from the flattened one
            call flattened_upper_triangle_index_to_element_indexes(itri, nstates, nstart, nstart, ie1, ie2)
            
            ! Avoid pairs with different degeneracy. They are part of different irreps,
            ! and the self-energy should be 0. The matter is discussed in detail in
            ! the subroutine add_q_contrib_to_offdiagonal_selfenergy_exchange_at_ik 
            ! within this file.
            if (band_degeneracy(ie1,ikp) == band_degeneracy(ie2,ikp)) then

                ! Loop over the intermediate band index 'm' in the formula
                do ie3 = mstart, mend
    
                  ! --- ENERGY CALCULATION (epsilon_{m, k-q}) ---
                  if (ie3 <= nstse) then
                    ! Standard valence or conduction band energy
                    energy_ie3 = evalfv(ie3, jkp)
                  else
                    ! Handling the core states
                    icg = ie3 - nstse
                    is  = corind(icg, 1)
                    ia  = corind(icg, 2)
                    ic  = corind(icg, 6)
                    ias = idxas(ia, is)
                    energy_ie3 = evalcr(ic, ias) - efermi
                  end if
    
                  
                  ! Complex energy variable
                  !   z_energy = ε_{m,k-q} - iω
                  z_energy = cmplx(energy_ie3, -freq_selfc%freqs(iom), dp)
    
                  ! Computes: [M^i_{nm}]* * W^c_{ij} * M^j_{lm}
                  ! wkq includes the normalization of the integration weight
                  mwm_offdiagonal = wkq * dot_product(minm(:,ie1,ie3), wm(:,ie2,ie3))
    
                  ! --- NUMERICAL CONVOLUTION (The Integral over omega') ---
                  ! Loop over the frequency grid 'iom2' which acts as omega'
                  do iom2 = 1, nomeg
    
                    if (iom2 == iom) then
                      ! --- RESIDUE TERM (Second line of the formula) ---
                      ! Adds the +1/2 * sgn(epsilon) * W(i*omega) term.
                      ! This occurs when the integration reaches the pole/specific frequency limit.
                      sigma_correlation_offdiagonal(itri,iom,ikp) = &
                          sigma_correlation_offdiagonal(itri,iom,ikp) + &
                          mwm_offdiagonal * sign(0.5_dp, energy_ie3)
                    else
                      ! We process only one source frequency (iom) at a time and therefore
                      ! only have access to the quantity
                      !
                      !     MWM(iom) = M† W(iω_iom) M .
                      !
                      ! The convolution contains the difference
                      !
                      !     W(iω') - W(iω),
                      !
                      ! so each stored MWM(iom) contributes to two different self-energy
                      ! values:
                      !
                      !   (1) a negative contribution to Σ(iω), corresponding to the
                      !       -W(iω) part of the integrand, and
                      !
                      !   (2) a positive contribution to Σ(iω'), corresponding to the
                      !       +W(iω') part when ω' = ω_iom.
                      !
                      ! Rather than storing MWM for all frequencies simultaneously, we
                      ! immediately accumulate both pieces:
                      !
                      !   Σ(iω)  -= K(ω,ω')  · MWM(iω)
                      !   Σ(iω') += K(ω',ω) · MWM(iω)
                      !
                      ! where K is the energy-dependent convolution kernel. After all
                      ! source frequencies have been processed, this exactly reconstructs
                      ! the term proportional to
                      !
                      !     W(iω') - W(iω),
                      !
                      ! while requiring only one frequency slice of MWM to be stored.
    
    
                      ! Complex energy variable
                      !   z_energy_p = ε_{m,k-q} - iω'
                      z_energy_p = cmplx(energy_ie3, -freq_selfc%freqs(iom2), dp)
    
                      ! Kernel evaluated for the current output frequency Σ(iω), using
                      ! ω' = freq(iom2) as the quadrature point.
                      energy_factor   = freq%womeg(iom2) / (freq%freqs(iom2)*freq%freqs(iom2) + z_energy*z_energy)
    
                      ! Companion kernel for the contribution deposited directly into
                      ! Σ(iω'), using the current frequency ω = freq(iom) as the
                      ! quadrature point.
                      energy_factor_p = freq%womeg(iom)  / (freq%freqs(iom)*freq%freqs(iom)  + z_energy_p*z_energy_p)
        
                      ! Contributions
                      convolution_contribution   = energy_factor   * z_energy   * mwm_offdiagonal / pi
                      convolution_contribution_p = energy_factor_p * z_energy_p * mwm_offdiagonal / pi
        
                      !   Σ(iω)  -= kernel(ω') * z(ω)  * MWM / π
                      sigma_correlation_offdiagonal(itri,iom,ikp) = &
                          sigma_correlation_offdiagonal(itri,iom,ikp) - convolution_contribution
        
                      !   Σ(iω') += kernel(ω) * z(ω') * MWM / π
                      sigma_correlation_offdiagonal(itri,iom2,ikp) = &
                          sigma_correlation_offdiagonal(itri,iom2,ikp) + convolution_contribution_p
                          
                    end if
    
                  end do ! iom2 (omega')
                end do ! ie3 (m)
            end if

        end do     ! itri
        !$omp end teams distribute parallel do
        OMP_OFFLOAD end target

        OMP_OFFLOAD target exit data map(delete: freq, freq_selfc)

    end subroutine add_q_omega_contrib_to_offdiagonal_selfenergy_correl_at_ik

    !> Read the offdiagonal part of correlation self-energy from files
    subroutine read_offdiagonal_selfenergy_correlation(kpt_indexes, file_format)
        use gw_io, only: build_file_name, read_bounds_from_file, read_from_file

        !> List of k-point indexes
        integer(i32), intent(in) :: kpt_indexes(:)
        !> Format of the file where to print. It can be e.g. 'text' or 'binary'
        character(len=*), intent(in) :: file_format

        integer(i32) :: i, lbounds(2), ubounds(2), nkpoints
        character(len=str_512) :: file_name
        nkpoints = size(kpt_indexes)

        call build_file_name( correlation_offdiagonal_selfenergy_rootname, kpt_indexes(1), file_name )
        call read_bounds_from_file( file_name, file_format, lbounds, ubounds )
        allocate( sigma_correlation_offdiagonal(lbounds(1):ubounds(1), lbounds(2):ubounds(2), kpt_indexes(1):kpt_indexes(nkpoints)))
        do i = 1, size( kpt_indexes )
            call build_file_name( correlation_offdiagonal_selfenergy_rootname, kpt_indexes(i), file_name )
            call read_from_file( file_name, sigma_correlation_offdiagonal(:,:,kpt_indexes(i)), lbounds, file_format )
        end do
    end subroutine read_offdiagonal_selfenergy_correlation

    !> Read the offdiagonal part of exchange self-energy from files
    subroutine read_offdiagonal_selfenergy_exchange(kpt_indexes, file_format)
        use gw_io, only: build_file_name, read_bounds_from_file, read_from_file

        !> List of k-point indexes
        integer(i32), intent(in) :: kpt_indexes(:)
        !> Format of the file where to print. It can be e.g. 'text' or 'binary'
        character(len=*), intent(in) :: file_format

        integer(i32) :: i, lbounds(1), ubounds(1), nkpoints 
        character(len=str_512) :: file_name

        nkpoints = size(kpt_indexes)

        call build_file_name( exchange_offdiagonal_selfenergy_rootname, kpt_indexes(1), file_name )
        call read_bounds_from_file( file_name, file_format, lbounds, ubounds )
        allocate( sigma_exchange_offdiagonal(lbounds(1):ubounds(1), kpt_indexes(1):kpt_indexes(nkpoints)))
        do i = 1, size( kpt_indexes )
            call build_file_name( exchange_offdiagonal_selfenergy_rootname, kpt_indexes(i), file_name )
            call read_from_file( file_name, sigma_exchange_offdiagonal(:, kpt_indexes(i)), lbounds(1), file_format )
        end do
    end subroutine read_offdiagonal_selfenergy_exchange

    !> Computes the off-diagonal GW self-energy matrix elements
    !>
    !>   Σ_{nl}(ε) = Re [Σ^x_{nl}
    !>             + 1/2 [ Σ^c_{nl}(ε_n) + Σ^c_{nl}(ε_l) ]]
    !>
    !> where Σ^c_{nl}(ω) is obtained by:
    !>   1) Analytic continuation from imaginary to real frequencies
    !>      using Padé approximants
    !>   2) Interpolation on the real axis using a 4th-order polynomial
    !>
    !> and adds them to the optimized xc non-local potential.
    !>
    !>
    !> Only band pairs belonging to the same irreducible representation
    !> (identified via band degeneracy) are considered.
    !> Diagonal terms are excluded.
    subroutine add_offdiag_selfenergy_at_energies_to_optimized_xc_potential(kpt_indexes, vxc_opt)

        use modgw,               only: kset
        use mod_selfenergy,      only: freq_selfc, generate_frequency_grid_for_correlation_self_energy
        use mod_bands,           only: evalfv
        use mod_frequency,       only: frequency, generate_freqgrid, delete_freqgrid
        use modinput,            only: input
        use mod_pade,            only: pade_approximant
        use mod_gw_degeneracies, only: band_degeneracy
        use modmpi,              only: terminate_if_false
        use math_utils,          only: flattened_upper_triangle_index_to_element_indexes, &
                                       number_of_upper_triangle_elements
                                       

    
        !> List of k-point indices where Σ is evaluated
        integer(i32), intent(in) :: kpt_indexes(:)
        !> Optimized potential
        complex(dp), allocatable, intent(inout) :: vxc_opt(:,:,:)
    
        integer(i32)       :: i, n, l, omega_idx
        integer(i32)       :: ni, nf, li, lf, num_l, num_n
        integer(long_int)  :: ntri, itri, offset
        real(dp)           :: omega, x
    
        complex(dp), allocatable :: sigmac_ac(:)
        complex(dp) :: dummy, sigma_n, sigma_l
    
        type(frequency) :: real_freq_selfc
    
        !------------------------------------------------------------------
        ! Generate real-frequency grid for analytic continuation and
        ! restore the self-energy grid.
        !------------------------------------------------------------------
        call generate_freqgrid(real_freq_selfc, &
                               input%gw%selfenergy%wgrid%type, &
                               'refreq', &
                               input%gw%selfenergy%wgrid%size, &
                               input%gw%selfenergy%wgrid%wmin, &
                               input%gw%selfenergy%wgrid%wmax)
        call delete_freqgrid(freq_selfc)
        call generate_frequency_grid_for_correlation_self_energy( input%gw )
    
        ni = lbound(vxc_opt, dim=1)
        nf = ubound(vxc_opt, dim=1)
        li = lbound(vxc_opt, dim=2)
        lf = ubound(vxc_opt, dim=2)

        !> Compute the number of upper triangle elements one have
        !> We only require the top half
        num_n   = nf - ni + 1_i32
        num_l   = lf - li + 1_i32
        call terminate_if_false( num_n == num_l, "Error(add_offdiagonal_selfenergy_at_given_energies_to_optimized_xc_potential): vxc is not a square matrix")
        ntri = number_of_upper_triangle_elements(num_n)

        !$omp parallel default(none) &
        !$omp shared(kpt_indexes, ni, li, band_degeneracy, sigma_correlation_offdiagonal) &
        !$omp shared(sigma_exchange_offdiagonal, freq_selfc, real_freq_selfc, evalfv, vxc_opt, ntri, num_n) &
        !$omp private(i, l, n, omega_idx, omega, sigmac_ac, sigma_n, sigma_l, dummy, itri, offset)
        allocate(sigmac_ac(real_freq_selfc%nomeg))
        !$omp do collapse(2)
        do i = 1, size(kpt_indexes)
            do itri = 1, ntri

                ! Compute the global indexes from the flattened one
                call flattened_upper_triangle_index_to_element_indexes(itri, num_n, ni, li, n, l)

                ! Skip bands belonging to different irreps
                if ( band_degeneracy(l, kpt_indexes(i)) /= &
                     band_degeneracy(n, kpt_indexes(i))) then
                     vxc_opt(n,l,kpt_indexes(i)) = zzero
                     vxc_opt(l,n,kpt_indexes(i)) = zzero
                else
                    sigmac_ac = zzero

                    !------------------------------------------------------
                    ! Analytic continuation to the real axis via Padé
                    !------------------------------------------------------
                    do omega_idx = 1, real_freq_selfc%nomeg
                        omega = real_freq_selfc%freqs(omega_idx)
    
                        if (omega < real_zero) then
                            call pade_approximant( &
                                freq_selfc%nomeg, &
                                cmplx(real_zero, -freq_selfc%freqs, dp), &
                                conjg(sigma_correlation_offdiagonal(itri,:,kpt_indexes(i))), &
                                cmplx(omega, real_zero, dp), &
                                sigmac_ac(omega_idx), &
                                dummy )
                        else
                            call pade_approximant( &
                                freq_selfc%nomeg, &
                                cmplx(real_zero, freq_selfc%freqs, dp), &
                                sigma_correlation_offdiagonal(itri,:,kpt_indexes(i)), &
                                cmplx(omega, real_zero, dp), &
                                sigmac_ac(omega_idx), &
                                dummy )
                        end if
                    end do
    
                    ! ------------------------------------------------------
                    ! Interpolate Σ^c(ω) at ε_n and ε_l (4th order)
                    ! ------------------------------------------------------
                    call get_selfc(real_freq_selfc%nomeg, real_freq_selfc%freqs, &
                                   sigmac_ac, evalfv(n, kpt_indexes(i)), &
                                   sigma_n, dummy)
    
                    call get_selfc(real_freq_selfc%nomeg, real_freq_selfc%freqs, &
                                   sigmac_ac, evalfv(l, kpt_indexes(i)), &
                                   sigma_l, dummy)
    
                    !------------------------------------------------------
                    ! Final off-diagonal self-energy contribution
                    ! to the optimized potential (real, Hermitian --- i.e.
                    ! a symmetric matrix)
                    !------------------------------------------------------
                    vxc_opt(n,l,kpt_indexes(i)) = &
                        real(0.5_dp * (sigma_l + sigma_n) &
                               + sigma_exchange_offdiagonal(itri,kpt_indexes(i)), kind=dp)
                    vxc_opt(l,n,kpt_indexes(i)) = vxc_opt(n,l,kpt_indexes(i))
                end if
            
            end do ! end upper tridiagonal elements loop
        end do ! end k loop
        !$omp end do
        deallocate(sigmac_ac)
        !$omp end parallel
    
        call delete_freqgrid(real_freq_selfc)
    
    end subroutine add_offdiag_selfenergy_at_energies_to_optimized_xc_potential

end module mod_offdiagonal_selfenergy
