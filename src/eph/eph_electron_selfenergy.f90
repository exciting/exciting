!> Module to compute electron-phonon contribution to electron self-energy
!> \(\Sigma_{nn'}({\bf k},\omega,T)\).
module eph_electron_selfenergy
  use eph_variables

  use precision, only: dp
#include "asserts.fpp"
  use modmpi

  implicit none
  private

  !> name for binary file to save electron self-energy
  character(*), parameter, public :: eph_else_filename = "EPH_ELSE.OUT"

  public :: eph_else_gen_fan_migdal_smearing, eph_else_gen_fan_migdal_aux, eph_else_gen_fan_migdal_from_aux, &
            eph_else_gen_debye_waller_ahc, &
            eph_else_hilo_setup_interpolation, eph_else_hilo_interpolate, &
            eph_else_gen_from_file, &
            eph_else_gen_specfun, eph_else_resample, &
            eph_else_solve_qp_equation, &
            eph_else_set_default_frequency_grid, eph_else_setting_string

contains

  !================================================================================ 
  ! FAN-MIGDAL SELF-ENERGY
  !
  !> Calculate auxiliary Fan-Migdal self-energy \(\Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\) using 
  !> tetrahedron integration.
  !>
  !> The auxiliary Fan-Migdal self-energy is given by
  !> \[
  !>    \Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T) = -\pi\, \sum\limits_{m,\nu} \int \frac{{\rm d}{\bf q}}{V_{\rm BZ}}
  !>    g_{mn,\nu}^\ast({\bf k},{\bf q})\, g_{mn',\nu}({\bf k},{\bf q})
  !>    \left[ \left( 1 - f_{m{\bf k}+{\bf q}}(T) + n_{\nu{\bf q}}(T) \right) \delta\left( \omega - \epsilon_{m{\bf k}+{\bf q}} - \omega_{\nu{\bf q}} \right)
  !>         + \left(     f_{m{\bf k}+{\bf q}}(T) + n_{\nu{\bf q}}(T) \right) \delta\left( \omega - \epsilon_{m{\bf k}+{\bf q}} + \omega_{\nu{\bf q}} \right) \right] \;,
  !> \]
  !> where the BZ integral is evaluated using tetrahedron integration.
  !>
  !> See [[eph_else_gen_fan_migdal_from_aux(subroutine)]] to see how to obtain the true Fan-Migdal self-energy from the auxiliary one.
  !>
  !> This subroutine does not include any internal MPI parallelism. But the sums over states \(m\) and phonon modes \(\nu\) can be 
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_gen_electron_selfenergy(subroutine)]].
  subroutine eph_else_gen_fan_migdal_aux( freqs, temps, el_energy_k, el_energy_kq, ph_energy_q, ephmat, tset, selfen_fm, &
      diagonal_only )
    use constants, only: zzero, pi
    use occupation_functions, only: fermi_dirac, bose_einstein
    use mod_opt_tetra, only: t_set, opt_tetra_int_delta
    use math_utils, only: get_degeneracies
    !> frequencies \(\omega\)
    real(dp), intent(in) :: freqs(:)
    !> temperatures \(\T\) in Kelvin
    real(dp), intent(in) :: temps(:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: el_energy_k(:)
    !> electron energies \(\epsilon_{m{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: el_energy_kq(:, :)
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy_q(:, :)
    !> EPH matrix elements \(g_{mn,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: ephmat(:, :, :, :)
    !> set of tetrahedra for integration
    type(t_set), intent(in) :: tset
    !> auxiliary Fan-Migdal self-energy \(\Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\)
    complex(dp), intent(out) :: selfen_fm(:, :, :)
    !> compute only band-diagonal elements \(\Sigma^{\rm FM}_{nn}({\bf k},\omega,T)\) only (default: `.true.`)
    logical, optional, intent(in) :: diagonal_only

    integer :: nwfk, nwfkq, nmode, nfreq, ntemp, nelem, nkq, iq, imode, itemp, ifreq, jst, i
    logical :: diag
    real(dp) :: fkq, nq

    integer, allocatable :: deg(:,:)
    real(dp), allocatable :: e(:,:)
    real(dp), allocatable :: g_squared(:), matr(:,:,:), intr(:,:,:)

    ! set defaults
    diag = .true.
    if (present(diagonal_only)) diag = diagonal_only

    ! set matrix sizes
    nwfk = size( ephmat, dim=2 )
    nwfkq = size( ephmat, dim=1 )
    nmode = size( ephmat, dim=3 )
    nkq = size( ephmat, dim=4 )
    nfreq = size( freqs )
    ntemp = size( temps )
    nelem = nwfk
    if (.not. diag) nelem = (nwfk * (nwfk + 1)) / 2 ! only upper triangle, Sigma^aux is Hermitian

    ! check input
    call terminate_if_false( diag, '(eph_else_gen_fan_migdal_aux): &
      Band-offdiagonal elments of self-energy not yet implemented.' )

    CALL_ASSERT( size( el_energy_k ) == nwfk, '1st dimension of `el_energy_k` must equal 2nd dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=1 ) == nwfkq, '1st dimension of `el_energy_kq` must equal 1st dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=2 ) == nkq, '2nd dimension of `el_energy_kq` must equal 4th dimension of `ephmat`.' )
    CALL_ASSERT( size( ph_energy_q, dim=1 ) == nmode, '1st dimension of `ph_energy_q` must equal 3rd dimension of `ephmat`.' )
    CALL_ASSERT( size( ph_energy_q, dim=2 ) == nkq, '2nd dimension of `ph_energy_q` must equal 4th dimension of `ephmat`.' )
    CALL_ASSERT( size( selfen_fm, dim=1 ) >= nfreq, '1st dimension of `selfen_fm` must be at least equal to length of array `freqs`.' )
    CALL_ASSERT( size( selfen_fm, dim=3 ) >= ntemp, '3rd dimension of `selfen_fm` must be at least equal to length of array `temps`.' )
    CALL_ASSERT( size( selfen_fm, dim=2 ) >= nelem, '2nd dimension of `selfen_fm` too small.' )

    selfen_fm(:nfreq, :nelem, :ntemp) = zzero

    deg = get_degeneracies( el_energy_k, eph_el_degtol )

    ! diagonal elements only (real integral)
    ! Note: For the band-diagonal elements, the auxiliary self-energy is equal to the imaginary part of the
    !       true self-energy. Hence, all integrals are real.
    if (diag) then
      allocate( e(nwfkq, nkq), g_squared(nelem), matr(nelem*ntemp, nwfkq, nkq), intr(nfreq, nelem*ntemp, nwfkq) )
      do imode = 1, nmode
        ! term with - \omega_{\nu q}
        !$omp parallel default( shared ) private( g_squared, itemp, fkq, nq )
        !$omp do collapse(2)
        do iq = 1, nkq
          do jst = 1, nwfkq
            g_squared = ephmat(jst, :, imode, iq)%re**2 + ephmat(jst, :, imode, iq)%im**2
            ! average over degenerate states
            do i = 1, size(deg, 2)
              g_squared(deg(1, i):deg(2, i)) = sum( g_squared(deg(1, i):deg(2, i)) ) / deg(3, i) 
            end do
            ! Note: delta(omega - (epsilon_{n k+q} + omega_{nu q})), hence '+'
            e(jst, iq) = el_energy_kq(jst, iq) + ph_energy_q(imode, iq)
            do itemp = 1, ntemp
              fkq = fermi_dirac( el_energy_kq(jst, iq), temps(itemp) )
              nq  = bose_einstein( ph_energy_q(imode, iq), temps(itemp) )
              matr(((itemp-1)*nelem+1):itemp*nelem, jst, iq) = (1.0_dp - fkq + nq) * g_squared
            end do
          end do
        end do
        !$omp end do
        !$omp end parallel        
        call opt_tetra_int_delta( tset, nkq, nwfkq, e, nfreq, freqs, &
          nelem*ntemp, intr, fun=matr )
        do jst = 1, nwfkq
          selfen_fm(1:nfreq, 1:nelem, 1:ntemp) = selfen_fm(1:nfreq, 1:nelem, 1:ntemp) &
            - cmplx( pi * reshape( intr(:, :, jst), [nfreq, nelem, ntemp] ), 0, dp )
        end do
        ! term with + \omega_{\nu q}
        !$omp parallel default( shared ) private( g_squared, itemp, fkq, nq )
        !$omp do collapse(2)
        do iq = 1, nkq
          do jst = 1, nwfkq
            g_squared = ephmat(jst, :, imode, iq)%re**2 + ephmat(jst, :, imode, iq)%im**2
            ! average over degenerate states
            do i = 1, size(deg, 2)
              g_squared(deg(1, i):deg(2, i)) = sum( g_squared(deg(1, i):deg(2, i)) ) / deg(3, i) 
            end do
            ! Note: delta(omega - (epsilon_{n k+q} - omega_{nu q})), hence '-'
            e(jst, iq) = el_energy_kq(jst, iq) - ph_energy_q(imode, iq)
            do itemp = 1, ntemp
              fkq = fermi_dirac( el_energy_kq(jst, iq), temps(itemp) )
              nq  = bose_einstein( ph_energy_q(imode, iq), temps(itemp) )
              matr(((itemp-1)*nelem+1):itemp*nelem, jst, iq) = (fkq + nq) * g_squared
            end do
          end do
        end do
        !$omp end do
        !$omp end parallel        
        call opt_tetra_int_delta( tset, nkq, nwfkq, e, nfreq, freqs, &
          nelem*ntemp, intr, fun=matr )
        do jst = 1, nwfkq
          selfen_fm(1:nfreq, 1:nelem, 1:ntemp) = selfen_fm(1:nfreq, 1:nelem, 1:ntemp) &
            - cmplx( pi * reshape( intr(:, :, jst), [nfreq, nelem, ntemp] ), 0, dp )
        end do
      end do
      deallocate( e, g_squared, matr, intr )
    end if
  end subroutine eph_else_gen_fan_migdal_aux

  !> Calculate Fan-Migdal self-energy \(\Sigma^{\rm FM}_{nn'}({\bf k},\omega,T)\) using 
  !> direct summation over BZ points \({\bf q}\) and imaginary smearing \(\eta\).
  !>
  !> The Fan-Migdal self-energy is given by
  !> \[
  !>    \Sigma^{\rm FM}_{nn'}({\bf k},\omega,T) = \sum\limits_{m,\nu} \int \frac{{\rm d}{\bf q}}{V_{\rm BZ}}
  !>    g_{mn,\nu}^\ast({\bf k},{\bf q})\, g_{mn',\nu}({\bf k},{\bf q})
  !>    \left[ \frac{1 - f_{m{\bf k}+{\bf q}}(T) + n_{\nu{\bf q}}(T)}{\omega - \epsilon_{m{\bf k}+{\bf q}} - \omega_{\nu{\bf q}} + {\rm i}\eta}
  !>         + \frac{    f_{m{\bf k}+{\bf q}}(T) + n_{\nu{\bf q}}(T)}{\omega - \epsilon_{m{\bf k}+{\bf q}} + \omega_{\nu{\bf q}} + {\rm i}\eta} \right] \;,
  !> \]
  !> where the BZ integral is approximated by a direct sum over a regular \({\bf q}\)-grid.
  !>
  !> This subroutine does not include any internal MPI parallelism. But the sums over states \(m\) and phonon modes \(\nu\) can be 
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_gen_electron_selfenergy(subroutine)]].
  subroutine eph_else_gen_fan_migdal_smearing( freqs, temps, el_energy_k, el_energy_kq, ph_energy_q, ephmat, qset, eta, selfen_fm, &
      diagonal_only )
    use constants, only: zzero
    use occupation_functions, only: fermi_dirac, bose_einstein
    use mod_kpointset, only: k_set
    use math_utils, only: get_degeneracies
    !> frequencies \(\omega\)
    real(dp), intent(in) :: freqs(:)
    !> temperatures \(\T\) in Kelvin
    real(dp), intent(in) :: temps(:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: el_energy_k(:)
    !> electron energies \(\epsilon_{m{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: el_energy_kq(:, :)
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy_q(:, :)
    !> EPH matrix elements \(g_{mn,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: ephmat(:, :, :, :)
    !> set of \({\bf q}\) vectors
    type(k_set), intent(in) :: qset
    !> imaginary broadening \(\eta\)
    real(dp), intent(in) :: eta
    !> Fan-Migdal self-energy \(\Sigma^{\rm FM}_{nn'}({\bf k},\omega,T)\)
    complex(dp), intent(out) :: selfen_fm(:, :, :)
    !> compute only band-diagonal elements \(\Sigma^{\rm FM}_{nn}({\bf k},\omega,T)\) only (default: `.true.`)
    logical, optional, intent(in) :: diagonal_only

    integer :: nwfk, nwfkq, nmode, nfreq, ntemp, nelem, iq, imode, itemp, ifreq, jst, i
    logical :: diag
    real(dp) :: fkq, nq

    integer, allocatable :: deg(:,:)
    real(dp), allocatable :: g_squared(:)

    ! set defaults
    diag = .true.
    if (present(diagonal_only)) diag = diagonal_only

    ! set matrix sizes
    nwfk = size( ephmat, dim=2 )
    nwfkq = size( ephmat, dim=1 )
    nmode = size( ephmat, dim=3 )
    nfreq = size( freqs )
    ntemp = size( temps )
    nelem = nwfk
    if( .not. diag ) nelem = nwfk**2

    ! check input
    call terminate_if_false( diag, '(eph_else_gen_fan_migdal_smearing): &
      Band-offdiagonal elements of self-energy not yet implemented.' )

    CALL_ASSERT( size( el_energy_k, dim=1 ) == nwfk,  '1st dimension of `el_energy_k` must equal 2nd dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=1 ) == nwfkq,  '1st dimension of `el_energy_kq` must equal 1st dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=2 ) == qset%nkpt,  '2nd dimension of `el_energy_kq` must equal number of points in `qset`.' )
    CALL_ASSERT( size( ph_energy_q, dim=1 ) == nmode,  '1st dimension of `ph_energy_q` must equal 3rd dimension of `ephmat`.' )
    CALL_ASSERT( size( ph_energy_q, dim=2 ) == qset%nkpt,  '2nd dimension of `ph_energy_q` must equal number of points `qset`.' )
    CALL_ASSERT( size( selfen_fm, dim=1 ) >= nfreq,  '1st dimension of `selfen_fm` must be at least equal to length of array `freqs`.' )
    CALL_ASSERT( size( selfen_fm, dim=3 ) >= ntemp,  '3rd dimension of `selfen_fm` must be at least equal to length of array `temps`.' )
    CALL_ASSERT( size( selfen_fm, dim=2 ) >= nelem,  '2nd dimension of `selfen_fm` too small.' )

    selfen_fm(:nfreq, :nelem, :ntemp) = zzero

    deg = get_degeneracies( el_energy_k, eph_el_degtol )

    ! diagonal elements only
    if (diag) then
      allocate( g_squared(nelem) )
      !$omp parallel default( shared ) private( g_squared, itemp, fkq, nq, ifreq ) reduction( +: selfen_fm )
      !$omp do collapse(3)
      do iq = 1, qset%nkpt
        do imode = 1, nmode
          do jst = 1, nwfkq
            g_squared = ephmat(jst, :, imode, iq)%re**2 + ephmat(jst, :, imode, iq)%im**2
            ! average over degenerate states
            do i = 1, size(deg, dim=2)
              g_squared(deg(1, i):deg(2, i)) = sum( g_squared(deg(1, i):deg(2, i)) ) / deg(3, i) 
            end do
            do itemp = 1, ntemp
              fkq = fermi_dirac( el_energy_kq(jst, iq), temps(itemp) )
              nq  = bose_einstein( ph_energy_q(imode, iq), temps(itemp) )
              do ifreq = 1, nfreq
                selfen_fm(ifreq, 1:nelem, itemp) = selfen_fm(ifreq, 1:nelem, itemp) + qset%wkpt(iq) * g_squared * &
                  ( (1.0_dp - fkq + nq) / cmplx( freqs(ifreq) - el_energy_kq(jst, iq) - ph_energy_q(imode, iq), eta, dp ) + &
                             (fkq + nq) / cmplx( freqs(ifreq) - el_energy_kq(jst, iq) + ph_energy_q(imode, iq), eta, dp ) )
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
      deallocate( g_squared )
    end if
  end subroutine eph_else_gen_fan_migdal_smearing

  !> Calculate true Fan-Migdal self-energy from auxiliary one as obtained by [[eph_else_gen_fan_migdal_aux(subroutine)]].
  !>
  !> The true Fan-Migdal self-energy \(\Sigma^{\rm FM}_{nn'}({\bf k},\omega,T)\) can be obtained from the 
  !> auxiliary self-energy \(\Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\) as
  !> \[ 
  !>    \Sigma^{\rm FM}_{nn'}({\bf k},\omega,T) = 
  !>    \left[ - \Re\, \tilde{\Sigma}^{\rm FM, aux}_{nn'}({\bf k},\omega,T) - \Im\, \Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T) \right]
  !>    + {\rm i} \left[ - \Im\, \tilde{\Sigma}^{\rm FM, aux}_{nn'}({\bf k},\omega,T) + \Re\, \Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T) \right] \;, 
  !>  \]
  !> where \(\tilde{\Sigma}^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\) is the Hilbert transform of \(\Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\).
  !>
  !> @note
  !> The Hilbert transform is computed using splines. Large oscillations in the spline fit can be suppressed by smoothing
  !> the input self-energy. Smooting can be achieved by a convolution with a Lorentzian of width `eta`. This is supposed
  !> to result in the same broadening as using the same value for `eta` in [[eph_else_gen_fan_migdal_smearing(subroutine)]].
  !> @endnote
  subroutine eph_else_gen_fan_migdal_from_aux( freqs, selfen_fm, swidth )
    use convolution, only: smoothen
    !> frequency grid \(\omega\) (in increasing order)
    real(dp), intent(in) :: freqs(:)
    !> on input: auxiliary Fan-Migdal self-energy \(\Simag^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\);   
    !> on output: true Fan-Migdal self-energy \(\Simag^{\rm FM}_{nn'}({\bf k},\omega,T)\)
    complex(dp), intent(inout) :: selfen_fm(:)
    !> target imaginary broadening \(\eta\) of Fan-Migdal self-energy (default: `0.0`)
    real(dp), optional, intent(in) :: swidth
    
    ! Note: For a more accurate result of the Hilbert transform, the frequency grid is extended to
    !       infinity on both sides.
    integer, parameter :: next = 30 ! number of extension frequency points

    integer :: nfreq, i
    real(dp) :: eta, df

    real(dp), allocatable :: freqs_ext(:), selfen_ext(:, :)

    ! set defaults
    eta = 0.0_dp
    if (present(swidth)) eta = swidth

    ! set matrix sizes
    nfreq = size( freqs )

    ! check input
    CALL_ASSERT( size( selfen_fm, dim=1 ) >= nfreq,  '1st dimension of `selfen_fm` must not be smaller than number of frequencies.' )
    CALL_ASSERT( nfreq >= 2,  'At least two frequency points are required.' )

    allocate( freqs_ext(nfreq+2*next) )
    allocate( selfen_ext(nfreq+2*next, 4), source=0.0_dp )

    ! generate extended frequency grid
    freqs_ext(next+1:next+nfreq) = freqs
    df = max( 1e-6_dp, freqs(2) - freqs(1) )
    do i = 0, next-1
      freqs_ext(next-i) = freqs(1) - 2**i * df
    end do
    df = max( 1e-6_dp, freqs(nfreq) - freqs(nfreq-1) )
    do i = 0, next-1
      freqs_ext(nfreq+next+i+1) = freqs(nfreq) + 2**i * df
    end do

    ! generate extended self-energy (zero outside given frequency range)
    selfen_ext(next+1:next+nfreq, 1) = selfen_fm(:nfreq)%re
    selfen_ext(next+1:next+nfreq, 2) = selfen_fm(:nfreq)%im

    ! apply Lorentzian broadening if necessary
    if (eta > epsilon(eta)) then
      call smoothen( freqs_ext, selfen_ext(:, 1), eta, kernel='lorentzian', mode='linear' )
      call smoothen( freqs_ext, selfen_ext(:, 2), eta, kernel='lorentzian', mode='linear' )
    end if

    ! apply Hilbert transform
    call hilbert_transform( nfreq+2*next, freqs_ext, selfen_ext(:, 1), selfen_ext(:, 3) )
    call hilbert_transform( nfreq+2*next, freqs_ext, selfen_ext(:, 2), selfen_ext(:, 4) )

    ! compose result
    selfen_fm(:nfreq) = cmplx( &
      - selfen_ext(next+1:next+nfreq, 3) - selfen_ext(next+1:next+nfreq, 2), &
      - selfen_ext(next+1:next+nfreq, 4) + selfen_ext(next+1:next+nfreq, 1), dp )

    deallocate( freqs_ext, selfen_ext )
  end subroutine eph_else_gen_fan_migdal_from_aux
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! DEBYE-WALLER SELF-ENERGY
  !
  !> Calculate band-diagonal Debye-Waller self-energy \(\Sigma^{\rm DW}_{nn}({\bf k},T)\) within
  !> the rigid-ion approximation (RIA) and Allen-Heine-Cardona (AHC) theory.
  !>
  !> The Debye-Waller self-energy is given by 
  !> \[ 
  !>    \Sigma^{\rm DW}_{nn}({\bf k},T) = \sum\limits_{\nu} \int \frac{{\rm d}{\bf q}}{V_{\rm BZ}} \sum\limits_{\kappa,\kappa',\alpha,\alpha'} 
  !>    \frac{e_{\kappa\alpha,\nu{\bf q}}^\ast}{\sqrt{2M_{\kappa}\, \omega_{\nu{\bf q}}}}\, 
  !>    g^{\rm DW}_{nn,\kappa\alpha,\kappa'\alpha'}({\bf k},{\bf q})\, 
  !>    \frac{e_{\kappa'\alpha',\nu{\bf q}}}{\sqrt{2M_{\kappa'}\, \omega_{\nu{\bf q}}}} 
  !>    \left(  n_{\nu{\bf q}}(T) + \frac{1}{2} \right) \;, 
  !> \]
  !> where the DW matrix elements are approximated as
  !> \[ 
  !>    g^{\rm DW}_{nn,\kappa\alpha,\kappa'\alpha'}({\bf k},{\bf q}) \approx - \delta_{\kappa,\kappa'}
  !>    \left( \sum\limits_{\kappa'',m} \frac{g_{mn,\kappa\alpha}^\ast({\bf k},{\bf \Gamma})\, g_{mn,\kappa''\alpha'}({\bf k},{\bf \Gamma})}{\epsilon_{n{\bf k}} - \epsilon_{m{\bf k}} + {\rm i}\eta} + c.c. \right)
  !>  \]
  !>
  !> This subroutine does not include any internal MPI parallelism. But the sums over states \(m\) and phonon modes \(\nu\) can be 
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_gen_electron_selfenergy(subroutine)]].
  subroutine eph_else_gen_debye_waller_ahc( temps, el_energy_n, el_energy_m, ph_energy_q, ph_evec, ephmat0, qset, eta, selfen_dw )
    use constants, only: zzero, zone
    use occupation_functions, only: bose_einstein
    use mod_kpointset, only: k_set
    use mod_atoms, only: natmtot, nspecies, natoms, idxas, spmass
    use modinput
    !> temperatures in Kelvin
    real(dp), intent(in) :: temps(:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: el_energy_n(:)
    !> electron energies \(\epsilon_{m{\bf k}}\)
    real(dp), intent(in) :: el_energy_m(:)
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy_q(:, :)
    !> phonon eigenvectors \(e_{\kappa\alpha,\nu}({\bf q})\)
    complex(dp), intent(in) :: ph_evec(:, :, :)
    !> EPH matrix elements \(g_{mn,\kappa\alpha}({\bf k},{\bf \Gamma})\) in atomic gauge
    complex(dp), intent(in) :: ephmat0(:, :, :)
    !> set of \({\bf q}\) vectors
    type(k_set), intent(in) :: qset
    !> imaginary broadening for degenerate states
    real(dp), intent(in) :: eta
    !> Debye-Waller self-energy \(\Sigma^{\rm DW}_{nn}({\bf k},T)\)
    complex(dp), intent(out) :: selfen_dw(:, :)
  
    integer :: nwfn, nwfm, nmode, ntemp, i, j, is, js, ia, ja, ias, jas, iq, imode, itemp, ist, jst
    real(dp) :: mi, de, nq

    ! variable names according to Eqs. (196-198)
    ! Giustino, F. Electron-phonon interactions from first principles. Rev. Mod. Phys. 89, 015003 (2017)
    complex(dp), allocatable :: t(:,:,:), ta(:,:,:), tsum(:,:), g(:,:,:)

    ! set matrix sizes
    nwfn = size( ephmat0, dim=2 )
    nwfm = size( ephmat0, dim=1 )
    nmode = size( ph_energy_q, dim=1 )
    ntemp = size( temps )

    ! check input
    CALL_ASSERT( size( ph_evec, dim=1 ) == 3*natmtot,  '1st dimension of `ph_evec` must equal 3 * number of atoms.' )
    CALL_ASSERT( size( ph_evec, dim=2 ) == nmode,  '2nd dimension of `ph_evec` must equal 1st dimension of `ph_energy_q`.' )
    CALL_ASSERT( size( ph_evec, dim=3 ) == qset%nkpt,  '3rd dimension of `ph_evec` must equal number of q-points.' )
    CALL_ASSERT( size( ph_energy_q, dim=2 ) == qset%nkpt,  '2nd dimension of `ph_energy_q` must equal number of q-points.' )
    CALL_ASSERT( size( el_energy_n, dim=1 ) == nwfn,  '1st dimension of `el_energy_n` must equal 2nd dimension of `ephmat0`.' )
    CALL_ASSERT( size( el_energy_m, dim=1 ) == nwfm,  '1st dimension of `el_energy_m` must equal 1st dimension of `ephmat0`.' )
    CALL_ASSERT( size( ephmat0, dim=3 ) == 3*natmtot,  '3rd dimension of `ephmat0` must equal 3 * number of atoms.' )
    CALL_ASSERT( size( selfen_dw, dim=2 ) >= ntemp,  '2nd dimension of `selfen_dw` must be at least equal to length of array `temps`.' )
    CALL_ASSERT( size( selfen_dw, dim=1 ) >= nwfn,  '1st dimension of `selfen_dw` too small.' )

    selfen_dw(1:nwfn, 1:ntemp) = zzero

    allocate( t(3*natmtot, 3*natmtot, qset%nkpt), ta(3, 3, natmtot), tsum(3*natmtot, 3*natmtot) )
    allocate( g(nwfm, nwfn, 3*natmtot) )

    do imode = 1, nmode
      t = zzero
      ! q-dependent part
      !$omp parallel default( shared ) private( ta, is, ia, ias, js, ja, jas, i, j, mi )
      !$omp do
      do iq = 1, qset%nkpt
        ta = zzero
        if (ph_energy_q(imode, iq) <= eph_ph_energy_zero) cycle
        do is = 1, nspecies
          if (spmass(is) == 0.0_dp) cycle
          mi = 1.0_dp / (2.0_dp * ph_energy_q(imode, iq) * spmass(is))
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            i = (ias - 1) * 3 + 1
            call zgerc( 3, 3, cmplx( mi, 0, dp ), ph_evec(i, imode, iq), 1, ph_evec(i, imode, iq), 1, ta(1, 1, ias), 3 )
          end do
        end do
        do js = 1, nspecies
          do ja = 1, natoms(js)
            jas = idxas(ja, js)
            j = (jas - 1) * 3
            do is = 1, nspecies
              do ia = 1, natoms(is)
                ias = idxas(ia, is)
                i = (ias - 1) * 3
                t(i+1:i+3, j+1:j+3, iq) = conjg( ta(:, :, ias) ) + ta(:, :, jas)
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      ! temperature dependent part
      do itemp = 1, ntemp
        tsum = zzero
        ! q-sum
        !$omp parallel default( shared ) private( nq ) reduction( +: tsum )
        !$omp do
        do iq = 1, qset%nkpt
          nq = bose_einstein( ph_energy_q(imode, iq), temps(itemp) )
          tsum = tsum + qset%wkpt(iq) * (nq + 0.5_dp) * t(:, :, iq)
        end do
        !$omp end do
        !$omp end parallel
        call zgemm( 'n', 't', nwfm*nwfn, 3*natmtot, 3*natmtot, zone, &
          ephmat0, nwfm*nwfn, tsum, 3*natmtot, zzero, g, nwfm*nwfn )
        g(:, :, 1) = sum( g * conjg( ephmat0 ), dim=3 )
        
        !$omp parallel default( shared ) private( de ) reduction( +: selfen_dw )
        !$omp do collapse(2)
        do ist = 1, nwfn
          do jst = 1, nwfm
            de = el_energy_n(ist) - el_energy_m(jst)
            selfen_dw(ist, itemp) = selfen_dw(ist, itemp) - g(jst, ist, 1) / cmplx( de, max(1e-6_dp, eta), dp )
          end do
        end do
        !$omp end do
        !$omp end parallel
      end do
    end do
  end subroutine eph_else_gen_debye_waller_ahc
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! HIGH AND LOW ENERGY SELF-ENERGY
  !
  !> Compute the high and low energy contribution to the electron self-energy in the static, adiabatic approximation.
  !>
  !> This subroutine does not include any internal MPI parallelism. But the sum over phonon modes \(\nu\) can be 
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_gen_electron_selfenergy(subroutine)]].
  subroutine eph_else_gen_hilo( temps, ph_energy_q, ephmat, qset, selfen_hilo )
    use constants, only: zzero, zone
    use occupation_functions, only: bose_einstein
    use mod_kpointset, only: k_set
    !> temperatures in Kelvin
    real(dp), intent(in) :: temps(:)
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy_q(:, :)
    !> hilo EPH matrix elements \(g^{\rm hilo}_{n,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: ephmat(:, :, :)
    !> set of \({\bf q}\) vectors
    type(k_set), intent(in) :: qset
    !> hilo self-energy \(\Sigma^{\rm hilo}_{nn}({\bf k},T)\)
    complex(dp), intent(out) :: selfen_hilo(:, :)
  
    integer :: nwfn, nmode, ntemp, iq, imode, itemp, ist
    real(dp) :: nq

    ! set matrix sizes
    nwfn = size( ephmat, dim=1 )
    nmode = size( ph_energy_q, dim=1 )
    ntemp = size( temps )

    ! check input
    CALL_ASSERT( size( ph_energy_q, dim=2 ) == qset%nkpt, '2nd dimension of `ph_energy_q` must equal number of q-points.' )
    CALL_ASSERT( size( ephmat, dim=2 ) == nmode, '2nd dimension of `ephmat` must equal 1st dimension of `ph_energy_q`.' )
    CALL_ASSERT( size( ephmat, dim=3 ) == qset%nkpt, '3rd dimension of `ephmat` must equal number of q-points.' )
    CALL_ASSERT( size( selfen_hilo, dim=2 ) >= ntemp, '2nd dimension of `selfen_hilo` must be at least equal to length of array `temps`.' )
    CALL_ASSERT( size( selfen_hilo, dim=1 ) >= nwfn, '1st dimension of `selfen_hilo` too small.' )

    selfen_hilo(1:nwfn, 1:ntemp) = zzero

    !$omp parallel default( shared ) private( nq ) reduction( +: selfen_hilo )
    !$omp do collapse(3)
    do imode = 1, nmode
      do itemp = 1, ntemp
        do iq = 1, qset%nkpt
          nq = bose_einstein( ph_energy_q(imode, iq), temps(itemp) )
          selfen_hilo(:, itemp) = selfen_hilo(:, itemp) + qset%wkpt(iq) * (nq + 0.5_dp) * ephmat(:, imode, iq)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine eph_else_gen_hilo

  !> Setup the interpolation of the high and low energy contribution to the electron self-energy.  
  !> The HILO self-energy can only be computed on the coarse electron \({\bf k}\)-grid. Its interpolation
  !> is based on the Wannier interpolation of the electron energies \(\epsilon_{n{\bf k}}\). 
  !> The HILO self-energy correction is assumed to be small compared to the bare electron energies.
  !> The renormalized electron energies on the coarse grid are interpolated to any arbitrary point
  !> and then the bare energy is subtracted from the renormalized one to obtain the HILO self-energy
  !> at the interpolation points.
  subroutine eph_else_hilo_setup_interpolation( temps, selfen_fm_R, selfen_dw_R )
    use eph_electrons, only: eph_el_energy_k, eph_el_evec_k, eph_el_mfi, eph_el_gen_Hk_wannier
    use eph_phonons, only: eph_ph_energy_q
    use eph_ephmat, only: eph_ephmat_hilo_read_coarse, eph_ephmat_gFMHILOkq_filename, eph_ephmat_gDWHILOkq_filename

    use block_data_file, only: block_data_file_type
    !> set of temperatures \(T\)
    real(dp), intent(in) :: temps(:)
    !> high and low energy contribution to Fan-Migdal self-energy in real-space Wannier gauge
    complex(dp), allocatable, intent(out) :: selfen_fm_R(:,:,:,:)
    !> high and low energy contribution to Debye-Waller self-energy in real-space Wannier gauge
    complex(dp), allocatable, intent(out) :: selfen_dw_R(:,:,:,:)
  
    integer :: ntemp, ik, iq, itemp, ik0, isym, ir
    type(block_data_file_type) :: gkq_file

    complex(dp), allocatable :: g(:,:,:), selfen_Hk(:,:), selfen_Wk(:,:,:,:)

    ntemp = size( temps )

    allocate( g(eph_fst:eph_lst, eph_nmode_tot, eph_qset_ph%nkpt) )
    allocate( selfen_Hk(eph_fst:eph_lst, ntemp) )
    allocate( selfen_Wk(eph_nwf_tot, eph_nwf_tot, eph_el_mfi%np, ntemp) )

    ! Fan-Migdal self-energy
    allocate( selfen_fm_R(eph_nwf_tot, eph_nwf_tot, eph_el_mfi%nr, ntemp) )
    gkq_file = block_data_file_type( eph_ephmat_gFMHILOkq_filename, [-1], cmplx( 0, 0, dp ) )
    call gkq_file%open( mpiglobal )
    do ik = 1, eph_el_mfi%np
      call findkptinset( eph_el_mfi%vpl(:, ik), eph_kset_el, isym, ik0 )
      ! read hilo matrix elements
      do iq = 1, eph_qset_ph%nkpt
        call eph_ephmat_hilo_read_coarse( eph_kset_el, eph_qset_ph, eph_el_mfi%vpl(:, ik), eph_qset_ph%vkl(:, iq), gkq_file, g(:, :, iq) )
      end do
      ! compute self-energy
      call eph_else_gen_hilo( temps, eph_ph_energy_q, g, eph_qset_ph, selfen_Hk )
      ! transform to Wannier gauge
      do itemp = 1, ntemp
        call eph_el_gen_Hk_wannier( eph_el_energy_k(:, ik0) + selfen_Hk(:, itemp)%re, eph_el_evec_k(:, :, ik0), selfen_Wk(:, :, ik, itemp) )
      end do
    end do
    ! transform to real space
    do itemp = 1, ntemp
      call eph_el_mfi%transform_p2R( [eph_nwf_tot, eph_nwf_tot], 1, selfen_Wk(:, :, :, itemp), eph_nwf_tot**2, 1, selfen_fm_R(:, :, :, itemp), eph_nwf_tot**2, 1 )
    end do
    call gkq_file%close( mpiglobal )

    ! Debye-Waller self-energy
    allocate( selfen_dw_R(eph_nwf_tot, eph_nwf_tot, eph_el_mfi%nr, ntemp) )
    gkq_file = block_data_file_type( eph_ephmat_gDWHILOkq_filename, [-1], cmplx( 0, 0, dp ) )
    call gkq_file%open( mpiglobal )
    do ik = 1, eph_el_mfi%np
      call findkptinset( eph_el_mfi%vpl(:, ik), eph_kset_el, isym, ik0 )
      ! read hilo matrix elements
      do iq = 1, eph_qset_ph%nkpt
        call eph_ephmat_hilo_read_coarse( eph_kset_el, eph_qset_ph, eph_el_mfi%vpl(:, ik), eph_qset_ph%vkl(:, iq), gkq_file, g(:, :, iq) )
      end do
      ! compute self-energy
      call eph_else_gen_hilo( temps, eph_ph_energy_q, g, eph_qset_ph, selfen_Hk )
      ! transform to Wannier gauge
      do itemp = 1, ntemp
        call eph_el_gen_Hk_wannier( eph_el_energy_k(:, ik0) + selfen_Hk(:, itemp)%re, eph_el_evec_k(:, :, ik0), selfen_Wk(:, :, ik, itemp) )
      end do
    end do
    ! transform to real space
    do itemp = 1, ntemp
      call eph_el_mfi%transform_p2R( [eph_nwf_tot, eph_nwf_tot], 1, selfen_Wk(:, :, :, itemp), eph_nwf_tot**2, 1, selfen_dw_R(:, :, :, itemp), eph_nwf_tot**2, 1 )
    end do
    call gkq_file%close( mpiglobal )

    deallocate( g, selfen_Hk, selfen_Wk )
  end subroutine eph_else_hilo_setup_interpolation

  !> Interpolate high and low energy contribution to the electron self-energy.  
  !> See [[eph_else_hilo_setup_interpolation(subroutine)]] for further information.
  subroutine eph_else_hilo_interpolate( vkl, el_energy_k, Umnk, selfen_R, selfen_k, &
      irange )
    use eph_electrons, only: eph_el_mfi, eph_el_mindist

    use constants, only: zzero, zone
    !> set of wave vectors \({\bf k}'\) in lattice coordinates
    real(dp), intent(in) :: vkl(:,:)
    !> electron energies \(\epsilon_{n{\bf k}}\) as obtained from [[eph_el_interpolate(subroutine)]]
    real(dp), intent(in) :: el_energy_k(:,:)
    !> Wannier transformation matrices \(U_{mn}({\bf k}')\) as obtained from [[eph_el_interpolate(subroutine)]]
    complex(dp), intent(in) :: Umnk(:,:,:)
    !> high and low energy self-energy in real-space Wannier gauge
    complex(dp), intent(in) :: selfen_R(:,:,:,:)
    !> high and low energy self-energy \(\Sigma^{rm hilo}_{nn}({\bf k},T)\)
    complex(dp), allocatable, intent(out) :: selfen_k(:,:,:)
    !> range of band index (default: all rows of \(U({\bf k})\))   
    !> Must match `irange` used in [[eph_el_interpolate(subroutine)]] to obtain `Umnk`. 
    integer, optional, intent(in) :: irange(2)

    integer :: irng(2), nk, nwfk, ntemp, ik, ist, itemp

    complex(dp), allocatable :: selfen_Wk(:,:,:), aux(:,:)

    complex(dp), external :: zdotu

    irng = [eph_fwf, eph_lwf]
    if (present(irange)) irng = irange

    ! set number of interpolation points
    nk = size( vkl, dim=2 )

    ! set matrix size
    nwfk = size( Umnk, dim=1 )
    ntemp = size( selfen_R, dim=4 )
    
    ! check input
    CALL_ASSERT( size( vkl, dim=1 ) == 3, '`vkl` must be a set of vectors of length 3.' )
    CALL_ASSERT( size( el_energy_k, dim=1 ) == nwfk, '1st dimension of `el_energy_k` must equal 1st dimension of `Umnk`.' )
    CALL_ASSERT( size( el_energy_k, dim=2 ) == nk, '2nd dimension of `el_energy_k` must equal number of k-vectors.' )
    CALL_ASSERT( size( Umnk, dim=3 ) == nk, '3rd dimension of `Umnk` must equal number of k-vectors.' )
    CALL_ASSERT( size( Umnk, dim=2 ) == eph_nwf_tot, '2nd dimension of `Umnk` must equal total number of Wannier functions.' )

    if (allocated(selfen_k)) deallocate( selfen_k )
    allocate( selfen_k(irng(1):irng(2), ntemp, nk) )

    allocate( selfen_Wk(eph_nwf_tot, eph_nwf_tot, nk), aux(eph_nwf_tot, nwfk) )

    do itemp = 1, ntemp
      ! interpolate self-energy
      call eph_el_mfi%transform_R2p( [eph_nwf_tot, eph_nwf_tot], 1, selfen_R(:, :, :, itemp), eph_nwf_tot**2, 1, selfen_Wk, eph_nwf_tot**2, 1, vkl, &
        minimal_distances=eph_el_mindist )
      ! transform self-energy to Hamiltonian gauge
      do ik = 1, nk
        call zgemm( 'n', 'c', eph_nwf_tot, nwfk, eph_nwf_tot, zone, &
          selfen_Wk(:, :, ik), eph_nwf_tot, &
          Umnk(:, :, ik), size(Umnk, dim=1), zzero, &
          aux, eph_nwf_tot )
        do ist = 1, nwfk
          selfen_k(irng(1)+ist-1, itemp, ik) = zdotu( eph_nwf_tot, Umnk(ist, 1, ik), size(Umnk, dim=1), aux(1, ist), 1 ) - el_energy_k(ist, ik)
        end do
      end do
    end do

    deallocate( selfen_Wk, aux )
  end subroutine eph_else_hilo_interpolate
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! SPECTRAL FUNCTION
  !
  !> Compute spectral function \(A_{nn'}({\bf k},\omega,T)\) from self-energy \(\Sigma_{nn'}({\bf k},\omega,T)\).
  !>
  !> The spectral function is defined as
  !> \[
  !>    A_n({\bf k},\omega,T) = - \frac{1}{\pi} \frac{\Im\, \Sigma_{nn}({\bf k},\omega,T)}
  !>    {\left(\omega - \epsilon_{n{\bf k}} - \Re\, \Sigma_{nn}({\bf k},\omega,T)\right)^2 + \left(\Im\, \Sigma_{nn}({\bf k},\omega,T)\right)^2} \;.
  !> \]
  subroutine eph_else_gen_specfun( freqs, selfen, sfun )
    use constants, only: pi
    !> shifted frequency grid \(\omega - \epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: freqs(:)
    !> self-energy \(\Sigma_{nn}({\bf k},\omega,T)\)
    complex(dp), intent(in) :: selfen(:)
    !> spectral function \(A_n({\bf k},\omega,T))
    real(dp), intent(out) :: sfun(:)
  
    real(dp), parameter :: eps_zero = 1e-12_dp
    real(dp), parameter :: piinv = -1.0_dp / pi

    integer :: n, i
    real(dp) :: eps
    
    n = size( freqs )

    ! check input
    CALL_ASSERT( size( selfen ) == n,  '`freqs` and `selfen` must be of equal size.' )
    CALL_ASSERT( size( sfun ) == n,  '`freqs` and `sfun` must be of equal size.' )

    ! infinitesimal to get peak even if Im(Sigma) is zero
    eps = sign( eps_zero, selfen(maxloc( abs( selfen%im ), dim=1 ))%im )

    sfun = piinv * (selfen%im + eps) / ((freqs - selfen%re)**2 + (selfen%im + eps)**2)
  end subroutine eph_else_gen_specfun
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! QUASI PARTICLE EQUATION
  !
  !> Solve quasi particle equation \(z = \epsilon^0 + \Sigma(z)\).
  !>
  !> Mode 0: on-the-mass-shell approximation   
  !>    \(\Sigma(z) \approx \Sigma(\epsilon^0)\)   
  !>    \(\epsilon^{\rm QP} = \epsilon^0 + Z\, \Sigma(\epsilon^0)\)   
  !>    \(Z = 1\)   
  !> Mode 1: real perturbative   
  !>    \(\Sigma(z) \approx \Sigma(\epsilon^0) + \Re \Sigma'(\epsilon^0)\, (z - \epsilon^0)\)   
  !>    \(\epsilon^{\rm QP} = \epsilon^0 + Z\, \Sigma(\epsilon^0)\)   
  !>    \(Z = \left( 1 - \Re \Sigma'(\epsilon^0) \right)^{-1}\)   
  !> Mode 2: complex perturbative   
  !>    \(\Sigma(z) \approx \Sigma(\epsilon^0) + \Sigma'(\epsilon^0)\, (z - \epsilon^0)\)   
  !>    \(\epsilon^{\rm QP} = \epsilon^0 + Z\, \Sigma(\epsilon^0)\)   
  !>    \(Z = \left( 1 - \Sigma'(\epsilon^0) \right)^{-1}\)   
  !> Mode 3: iterative
  !>    \(\Re \epsilon^{\rm QP} = \epsilon^0 + \Re \Sigma(\Re \epsilon^{\rm QP})\)   
  !>    \(\epsilon^{\rm QP} = \epsilon^0 + \Sigma(\Re \epsilon^{\rm QP})\)
  subroutine eph_else_solve_qp_equation( e0, freqs, selfen, eQP, selfen0, Z, mode )
    use math_utils, only: interp1d
    use savitzky_golay, only: savgol
    !> single-particle energy \(\epsilon^0\)
    real(dp), intent(in) :: e0
    !> real frequencies \(\omega\)
    real(dp), intent(in) :: freqs(:)
    !> self-energy \(\Sigma(\omega)\)
    complex(dp), intent(in) :: selfen(:)
    !> quasi-particle energy \(\epsilon^{\rm QP}\)
    complex(dp), intent(out) :: eQP
    !> self-energy at evaluation point (\(\epsilon^0\) for modes 0,1,2 and \(\Re \epsilon^{\rm QP}\) for mode 3)
    complex(dp), intent(out) :: selfen0
    !> quasi-particle strength \(Z\)
    complex(dp), intent(out) :: Z
    !> solver mode
    integer, intent(in) :: mode

    integer, parameter :: MAX_ITER = 100 ! maximum number of iterations
    
    integer :: nfreq, i
    real(dp) :: re0(2), im0(2)

    real(dp), allocatable :: re(:,:), im(:,:)

    nfreq = size(freqs)

    allocate( re(nfreq, 2), im(nfreq, 2) )

    re(:, 1) = selfen%re
    im(:, 1) = selfen%im

    ! compute smoothed derivative
    call savgol( freqs, re, degree=1, window=3, derivative=[0, 1] )
    call savgol( freqs, im, degree=1, window=3, derivative=[0, 1] )

    select case (mode)
      case (0, 1, 2)
        call interp1d( freqs, re(:, 1), [e0], re0(1:1), 'spline' )
        call interp1d( freqs, re(:, 2), [e0], re0(2:2), 'spline' )
        call interp1d( freqs, im(:, 1), [e0], im0(1:1), 'spline' )
        call interp1d( freqs, im(:, 2), [e0], im0(2:2), 'spline' )
        select case (mode)
          case (0)
            Z = cmplx( 1, 0, dp) 
          case (1)
            Z = cmplx( 1, 0, dp) / cmplx( 1.0_dp - re0(2), 0, dp ) 
          case (2)
            Z = cmplx( 1, 0, dp) / cmplx( 1.0_dp - re0(2), -im0(2), dp ) 
          case default
            Z = cmplx( 0, 0, dp )
        end select
      case (3)
        re0(1) = e0
        do i = 1, MAX_ITER
          call interp1d( freqs, re(:, 1), re0(1:1), re0(2:2), 'spline' )
          re0(2) = e0 + re0(2)
          if (abs( re0(2) - re0(1) ) < 2*epsilon(re0)) exit 
          re0(1) = re0(2)
        end do
        call interp1d( freqs, im(:, 1), re0(1:1), im0(1:1), 'spline' )
        re0(1) = re0(1) - e0
        Z = cmplx( 1, 0, dp ) 
      case default
        call terminate_if_false( .false., '(eph_else_solve_qp_equation) &
          Invalid solver mode.' )
    end select

    selfen0 = cmplx( re0(1), im0(1), dp )
    eQP = e0 + Z*selfen0
  end subroutine eph_else_solve_qp_equation
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! AUXILIARY PROCEDURES
  !
  !> Read self-energy from file and apply postprocessing if necessary.
  subroutine eph_else_gen_from_file( file, ik, freqs, temps, el_energy, selfen_fm, selfen_dw, &
      selfen_fm_hilo, selfen_dw_hilo, swidth )
    use eph_inout, only: eph_io_read_el_self_energy
    use block_data_file, only: block_data_file_type
    use convolution, only: smoothen
    !> binary file
    type(block_data_file_type), intent(inout) :: file
    !> index of \({\bf k}\)-point
    integer, intent(in) :: ik
    !> frequencies \(\omega\)
    real(dp), allocatable, intent(out) :: freqs(:)
    !> temperatures \(T\)
    real(dp), allocatable, intent(out) :: temps(:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), allocatable, intent(out) :: el_energy(:)
    !> Fan-Migdal self-energy \(\Sigma^{\rm FM}_n({\bf k},\omega,T)\)
    complex(dp), allocatable, intent(out) :: selfen_fm(:,:,:)
    !> Debye-Waller self-energy \(\Sigma^{\rm DW}_n({\bf k},T)\)
    real(dp), allocatable, intent(out) :: selfen_dw(:,:)
    !> high and low energy Fan-Migdal self-energy \(\Sigma^{\rm FM, hilo}_n({\bf k},T)\)
    real(dp), allocatable, optional, intent(out) :: selfen_fm_hilo(:,:)
    !> high and low energy Debye-Waller self-energy \(\Sigma^{\rm DW, hilo}_n({\bf k},T)\)
    real(dp), allocatable, optional, intent(out) :: selfen_dw_hilo(:,:)
    !> target imaginary broadening \(\eta\) of Fan-Migdal self-energy (default: `0.0`)
    real(dp), optional, intent(in) :: swidth
    
    integer :: itemp, ist, fst
    real(dp) :: eta, eta0
    character(len=:), allocatable :: integration

    real(dp), allocatable :: fun1(:), fun2(:), fm_hilo(:,:), dw_hilo(:,:)

    ! set defaults
    eta = 0.0_dp
    if (present(swidth)) eta = swidth

    ! read data from file
    call eph_io_read_el_self_energy( file, ik, fst, freqs, temps, el_energy, selfen_fm, selfen_dw, fm_hilo, dw_hilo, integration, eta0 )

    ! apply Lorentzian broadening if necessary
    eta = eta - eta0 
    if (eta > epsilon(eta) .and. (integration /= 'kramers-kronig')) then
      do itemp = 1, size(temps)
        do ist = lbound(el_energy, dim=1), ubound(el_energy, dim=1)
          fun1 = selfen_fm(:, ist, itemp)%re
          fun2 = selfen_fm(:, ist, itemp)%im
          call smoothen( freqs, fun1, eta, kernel='lorentzian', mode='linear' )
          call smoothen( freqs, fun2, eta, kernel='lorentzian', mode='linear' )
          selfen_fm(:, ist, itemp) = cmplx( fun1, fun2, dp )
        end do
      end do
    end if

    ! obtain true self-energy from auxiliary one
    if (integration == 'kramers-kronig') then
      do itemp = 1, size(temps)
        do ist = lbound(el_energy, dim=1), ubound(el_energy, dim=1)
          call eph_else_gen_fan_migdal_from_aux( freqs, selfen_fm(:, ist, itemp), swidth=eta )
        end do
      end do
    end if

    ! add high and low energy contribution
    do itemp = 1, size(temps)
      do ist = lbound(el_energy, dim=1), ubound(el_energy, dim=1)
        selfen_fm(:, ist, itemp) = selfen_fm(:, ist, itemp) + fm_hilo(ist, itemp) 
        selfen_dw(ist, itemp) = selfen_dw(ist, itemp) + dw_hilo(ist, itemp) 
      end do
    end do

    if (present(selfen_fm_hilo)) selfen_fm_hilo = fm_hilo
    if (present(selfen_dw_hilo)) selfen_dw_hilo = dw_hilo
  end subroutine eph_else_gen_from_file

  !> Set default parameters for frequency grid used for electron self-energy calculation.
  pure subroutine eph_else_set_default_frequency_grid( fgrid )
    use modinput, only: freq_grid_type
    !> frequency grid object
    type(freq_grid_type), intent(out) :: fgrid
  
    fgrid%type = 'density'          ! density based sampling
    fgrid%numpoints = 500           ! number of sampling points
    fgrid%range = [0.0_dp, 0.0_dp]  ! automatic range detection
    fgrid%padding = 0.5_dp          ! padding to add at both ends of range
    fgrid%lorentzwidth = 0.04_dp    ! width of Lorantzian density
  end subroutine eph_else_set_default_frequency_grid

  !> Get string representing settings used for self-energy calculation.
  function eph_else_setting_string( ngridbz, wfrange, nfreq, temps ) result( settings )
    !> Brillouin zone sampling used for integration
    integer, intent(in) :: ngridbz(3)
    !> Wannier function range
    integer, intent(in) :: wfrange(2)
    !> number of frequencies
    integer, intent(in) :: nfreq
    !> temperatures
    real(dp), intent(in) :: temps(:)
    !> setting string
    character(:), allocatable :: settings

    integer :: ntemp, tmin, tmax, tstep
    character(8) :: bz_digit, wf_digit, freq_digit, temp_digit
    character(128) :: string

    ntemp = size(temps)
    tmin = nint( temps(1) )
    tmax = nint( temps(ntemp) )
    tstep = 0
    if (ntemp > 1) tstep = nint( abs( temps(ntemp) - temps(1) ) / (ntemp - 1) )

    ! set number of digits
    write( bz_digit, '(i8.8)' ) max( 2, ceiling( log10( maxval( ngridbz ) + 1e-12_dp ) ) )
    write( wf_digit, '(i8.8)' ) max( 2, ceiling( log10( maxval( wfrange ) + 1e-12_dp ) ) )
    write( freq_digit, '(i8.8)' ) max( 3, ceiling( log10( nfreq + 1e-12_dp ) ) )
    write( temp_digit, '(i8.8)' ) max( 3, ceiling( log10( maxval( [tmin, tmax, tstep] ) + 1e-12_dp ) ) )

    write( string, '(&
      &"BZ",3(i'//bz_digit//'.'//bz_digit//',"_"),&
      &"WF",2(i'//wf_digit//'.'//wf_digit//',"_"),&
      &"F",i'//freq_digit//'.'//freq_digit//',"_",&
      &"T",2(i'//temp_digit//'.'//temp_digit//',"_"),i'//temp_digit//'.'//temp_digit//')' ) &
      ngridbz, wfrange, nfreq, tmin, tmax, tstep

    settings = trim( string )
  end function eph_else_setting_string

  !> Resample self-energy on new frequency grid.
  !>
  !> If `fgrid_out%type == 'density'`, the output sampling is such, that the self-energy is densely
  !> sampled around the features in the spectral function.
  subroutine eph_else_resample( nfreq_in, freqs_in, selfen_in, fgrid_out, freqs_out, selfen_out, &
      interpolation_method )
    use grid_utils, only: linspace, spacing_from_density
    use math_utils, only: interp1d
    use convolution, only: smoothen
    use modinput, only: freq_grid_type
    !> number frequency points in input data
    integer, intent(in) :: nfreq_in
    !> input frequency grid
    real(dp), intent(in) :: freqs_in(nfreq_in)
    !> input self-energy
    complex(dp), intent(in) :: selfen_in(nfreq_in)
    !> frequency grid for resampled data
    type(freq_grid_type), intent(in) :: fgrid_out
    !> resampled frequency grid
    real(dp), intent(out) :: freqs_out(fgrid_out%numpoints)
    !> resampled self-energy
    complex(dp), intent(out) :: selfen_out(fgrid_out%numpoints)
    !> interpolation mode (`'linear'` or `'spline'`)    
    !> default: `'spline'`
    character(len=*), optional, intent(in) :: interpolation_method
  
    ! weight of spectral function in point density
    ! 0 - uniform sampling
    ! 1 - no points, where spectral function is zero
    real(dp), parameter :: sfun_wgt = 0.9_dp

    real(dp), allocatable :: sfun_in(:), point_density(:), fun1(:), fun2(:)
    character(len=:), allocatable :: method

    method = 'spline'
    if (present(interpolation_method)) method = trim( adjustl( interpolation_method ) )

    allocate( fun2(fgrid_out%numpoints), sfun_in(nfreq_in) )

    ! generate output sampling
    fun1 = linspace( fgrid_out%range(1)-fgrid_out%padding, fgrid_out%range(2)+fgrid_out%padding, fgrid_out%numpoints )
    select case (fgrid_out%type)
      case ('uniform')
        freqs_out = fun1
      case ('density')
        ! compute spectral function on input grid
        call eph_else_gen_specfun( freqs_in, selfen_in, sfun_in )
        ! generate point density
        point_density = sfun_wgt * abs(sfun_in) + (1.0_dp - sfun_wgt) / (freqs_in(nfreq_in) - freqs_in(1))
        call smoothen( freqs_in, point_density, fgrid_out%lorentzwidth, kernel='lorentzian', mode='linear' )
        call interp1d( freqs_in, point_density, fun1, fun2, method='linear' )
        ! generate output sampling according to point density
        freqs_out = spacing_from_density( fun1, fun2, fgrid_out%numpoints )
    end select
    ! interpolate self-energy to output sampling
    select case (method)
      case ('spline')
        call interp1d( freqs_in, (selfen_in%re), freqs_out, fun1, method='spline' )
        call interp1d( freqs_in, (selfen_in%im), freqs_out, fun2, method='spline' )
      case ('linear')
        call interp1d( freqs_in, (selfen_in%re), freqs_out, fun1, method='linear' )
        call interp1d( freqs_in, (selfen_in%im), freqs_out, fun2, method='linear' )
      case default
        call terminate_if_false( .false., '(eph_else_resample) &
          Unsupported interpolation mode `'//method//'`.' )
    end select
    selfen_out = cmplx( fun1, fun2, dp )
  end subroutine eph_else_resample
  !-------------------------------------------------------------------------------- 

end module eph_electron_selfenergy
