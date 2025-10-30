!> Module to compute electron-phonon contribution to electron self-energy
!> \(\Sigma_{nn'}({\bf k},\omega,T)\).
module eph_electron_selfenergy
  use eph_variables

  use precision, only: dp
  use asserts, only: assert
  use modmpi

  implicit none
  private

  public :: eph_else_set_default_frequency_grid, eph_else_get_binary_file_name, &
            eph_else_gen_fan_migdal_smearing, eph_else_gen_fan_migdal_aux, eph_else_gen_fan_migdal_from_aux, &
            eph_else_gen_debye_waller_ahc, &
            eph_else_gen_specfun, eph_else_resample

contains

  !================================================================================ 
  ! FAN-MIGDAL SELF-ENERGY
  !
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
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_electron_selfenergy_gen(subroutine)]].
  subroutine eph_else_gen_fan_migdal_smearing( freqs, temps, el_energy, ph_energy, ephmat, qset, eta, selfen_fm, &
      diagonal_only )
    use constants, only: zzero
    use occupation_functions, only: fermi_dirac, bose_einstein
    use mod_kpointset, only: k_set
    !> frequencies \(\omega\)
    real(dp), intent(in) :: freqs(:)
    !> temperatures \(\T\) in Kelvin
    real(dp), intent(in) :: temps(:)
    !> electron energies \(\epsilon_{m{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: el_energy(:, :)
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy(:, :)
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

    integer :: nwfk, nwfkq, nmode, nfreq, ntemp, nelem, iq, imode, itemp, ifreq, jst
    logical :: diag
    real(dp) :: fkq, nq

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

    call assert( size( el_energy, dim=1 ) == nwfkq, &
      '1st dimension of `el_energy` must equal 1st dimension of `ephmat`.' )
    call assert( size( el_energy, dim=2 ) == qset%nkpt, &
      '2nd dimension of `el_energy` must equal number of points in `qset`.' )
    call assert( size( ph_energy, dim=1 ) == nmode, &
      '1st dimension of `ph_energy` must equal 3rd dimension of `ephmat`.' )
    call assert( size( ph_energy, dim=2 ) == qset%nkpt, &
      '2nd dimension of `ph_energy` must equal number of points `qset`.' )
    call assert( size( selfen_fm, dim=1 ) >= nfreq, &
      '1st dimension of `selfen_fm` must be at least equal to length of array `freqs`.' )
    call assert( size( selfen_fm, dim=3 ) >= ntemp, &
      '3rd dimension of `selfen_fm` must be at least equal to length of array `temps`.' )
    call assert( size( selfen_fm, dim=2 ) >= nelem, &
      '2nd dimension of `selfen_fm` too small.' )

    selfen_fm(:nfreq, :nelem, :ntemp) = zzero

    ! diagonal elements only
    if (diag) then
      allocate( g_squared(nelem) )
      !$omp parallel default( shared ) private( g_squared, itemp, fkq, nq, ifreq ) reduction( +: selfen_fm )
      !$omp do collapse(3)
      do iq = 1, qset%nkpt
        do imode = 1, nmode
          do jst = 1, nwfkq
            g_squared = ephmat(jst, :, imode, iq)%re**2 + ephmat(jst, :, imode, iq)%im**2
            do itemp = 1, ntemp
              fkq = fermi_dirac( el_energy(jst, iq), temps(itemp) )
              nq  = bose_einstein( ph_energy(imode, iq), temps(itemp) )
              do ifreq = 1, nfreq
                selfen_fm(ifreq, 1:nelem, itemp) = selfen_fm(ifreq, 1:nelem, itemp) + qset%wkpt(iq) * g_squared * &
                  ( (1.0_dp - fkq + nq) / cmplx( freqs(ifreq) - el_energy(jst, iq) - ph_energy(imode, iq), eta, dp ) + &
                             (fkq + nq) / cmplx( freqs(ifreq) - el_energy(jst, iq) + ph_energy(imode, iq), eta, dp ) )
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
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_electron_selfenergy_gen(subroutine)]].
  subroutine eph_else_gen_fan_migdal_aux( freqs, temps, el_energy, ph_energy, ephmat, tset, selfen_fm, &
      diagonal_only )
    use constants, only: zzero, pi
    use occupation_functions, only: fermi_dirac, bose_einstein
    use mod_opt_tetra, only: t_set, opt_tetra_int_delta
    !> frequencies \(\omega\)
    real(dp), intent(in) :: freqs(:)
    !> temperatures \(\T\) in Kelvin
    real(dp), intent(in) :: temps(:)
    !> electron energies \(\epsilon_{m{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: el_energy(:, :)
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy(:, :)
    !> EPH matrix elements \(g_{mn,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: ephmat(:, :, :, :)
    !> set of tetrahedra for integration
    type(t_set), intent(in) :: tset
    !> auxiliary Fan-Migdal self-energy \(\Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\)
    complex(dp), intent(out) :: selfen_fm(:, :, :)
    !> compute only band-diagonal elements \(\Sigma^{\rm FM}_{nn}({\bf k},\omega,T)\) only (default: `.true.`)
    logical, optional, intent(in) :: diagonal_only

    integer :: nwfk, nwfkq, nmode, nfreq, ntemp, nelem, nkq, iq, imode, itemp, ifreq, jst
    logical :: diag
    real(dp) :: fkq, nq

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

    call assert( size( el_energy, dim=1 ) == nwfkq, &
      '1st dimension of `el_energy` must equal 1st dimension of `ephmat`.' )
    call assert( size( el_energy, dim=2 ) == nkq, &
      '2nd dimension of `el_energy` must equal 4th dimension of `ephmat`.' )
    call assert( size( ph_energy, dim=1 ) == nmode, &
      '1st dimension of `ph_energy` must equal 3rd dimension of `ephmat`.' )
    call assert( size( ph_energy, dim=2 ) == nkq, &
      '2nd dimension of `ph_energy` must equal 4th dimension of `ephmat`.' )
    call assert( size( selfen_fm, dim=1 ) >= nfreq, &
      '1st dimension of `selfen_fm` must be at least equal to length of array `freqs`.' )
    call assert( size( selfen_fm, dim=3 ) >= ntemp, &
      '3rd dimension of `selfen_fm` must be at least equal to length of array `temps`.' )
    call assert( size( selfen_fm, dim=2 ) >= nelem, &
      '2nd dimension of `selfen_fm` too small.' )

    selfen_fm(:nfreq, :nelem, :ntemp) = zzero

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
            ! Note: delta(omega - (epsilon_{n k+q} + omega_{nu q})), hence '+'
            e(jst, iq) = el_energy(jst, iq) + ph_energy(imode, iq)
            do itemp = 1, ntemp
              fkq = fermi_dirac( el_energy(jst, iq), temps(itemp) )
              nq  = bose_einstein( ph_energy(imode, iq), temps(itemp) )
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
            ! Note: delta(omega - (epsilon_{n k+q} - omega_{nu q})), hence '-'
            e(jst, iq) = el_energy(jst, iq) - ph_energy(imode, iq)
            do itemp = 1, ntemp
              fkq = fermi_dirac( el_energy(jst, iq), temps(itemp) )
              nq  = bose_einstein( ph_energy(imode, iq), temps(itemp) )
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

  !> Calculate true Fan-Migdal self-energy from auxiliary one as obtained by [[eph_else_gen_fan_migdal_aux(subroutine)]].
  !>
  !> The true Fan-Migdal self-energy \(\Sigma^{\rm FM}_{nn'}({\bf k},\omega,T)\) can be obtained from the 
  !> auxiliary self-energy \(\Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\) as
  !> \[ 
  !>    \Sigma^{\rm FM}_{nn'}({\bf k},\omega,T) = 
  !>    \left[ \Re\, \tilde{\Sigma}^{\rm FM, aux}_{nn'}({\bf k},\omega,T) - \Im\, \Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T) \right]
  !>    + {\rm i} \left[ \Im\, \tilde{\Sigma}^{\rm FM, aux}_{nn'}({\bf k},\omega,T) + \Re\, \Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T) \right] \;, 
  !>  \]
  !> where \(\tilde{\Sigma}^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\) is the Hilbert transform of \(\Sigma^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\).
  !>
  !>
  !> The Hilbert transform is computed using splines. Large oscillations in the spline fit can be suppressed by smoothing
  !> the input self-energy. Smooting can be achieved by a convolution with a Lorentzian of width `eta`. This is supposed
  !> to result in the same broadening as using the same value for `eta` in [[eph_else_gen_fan_migdal_smearing(subroutine)]].
  !>
  subroutine eph_else_gen_fan_migdal_from_aux( freqs, selfen_fm )
    !> frequency grid \(\omega\) (in increasing order)
    real(dp), intent(in) :: freqs(:)
    !> on input: auxiliary Fan-Migdal self-energy \(\Simag^{\rm FM, aux}_{nn'}({\bf k},\omega,T)\);   
    !> on output: true Fan-Migdal self-energy \(\Simag^{\rm FM}_{nn'}({\bf k},\omega,T)\)
    complex(dp), intent(inout) :: selfen_fm(:)
    
    ! Note: For a more accurate result of the Hilbert transform, the frequency grid is extended to
    !       infinity on both sides.
    integer, parameter :: next = 30 ! number of extension frequency points

    integer :: nfreq, i
    real(dp) :: swidth, df

    real(dp), allocatable :: freqs_ext(:), selfen_ext(:, :)

    ! set matrix sizes
    nfreq = size( freqs )

    ! check input
    call assert( size( selfen_fm, dim=1 ) >= nfreq, &
      '1st dimension of `selfen_fm` must not be smaller than number of frequencies.' )
    call assert( nfreq >= 2, &
      'At least two frequency points are required.' )

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

    ! apply Hilbert transform
    call hilbert_transform( nfreq+2*next, freqs_ext, selfen_ext(:, 1), selfen_ext(:, 3) )
    call hilbert_transform( nfreq+2*next, freqs_ext, selfen_ext(:, 2), selfen_ext(:, 4) )

    ! compose result
    selfen_fm(:nfreq) = cmplx( &
      selfen_ext(next+1:next+nfreq, 3) - selfen_ext(next+1:next+nfreq, 2), &
      selfen_ext(next+1:next+nfreq, 4) + selfen_ext(next+1:next+nfreq, 1), dp )

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
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_electron_selfenergy_gen(subroutine)]].
  subroutine eph_else_gen_debye_waller_ahc( temps, el_energy_n, el_energy_m, ph_energy, ph_evec, ephmat0, qset, eta, selfen_dw )
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
    real(dp), intent(in) :: ph_energy(:, :)
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
    nmode = size( ph_energy, dim=1 )
    ntemp = size( temps )

    ! check input
    call assert( size( ph_evec, dim=1 ) == 3*natmtot, &
      '1st dimension of `ph_evec` must equal 3 * number of atoms.' )
    call assert( size( ph_evec, dim=2 ) == nmode, &
      '2nd dimension of `ph_evec` must equal 1st dimension of `ph_energy_q`.' )
    call assert( size( ph_evec, dim=3 ) == qset%nkpt, &
      '3rd dimension of `ph_evec` must equal number of q-points.' )
    call assert( size( ph_energy, dim=2 ) == qset%nkpt, &
      '2nd dimension of `ph_energy` must equal number of q-points.' )
    call assert( size( el_energy_n, dim=1 ) == nwfn, &
      '1st dimension of `el_energy_n` must equal 2nd dimension of `ephmat0`.' )
    call assert( size( el_energy_m, dim=1 ) == nwfm, &
      '1st dimension of `el_energy_m` must equal 1st dimension of `ephmat0`.' )
    call assert( size( ephmat0, dim=3 ) == 3*natmtot, &
      '3rd dimension of `ephmat0` must equal 3 * number of atoms.' )
    call assert( size( selfen_dw, dim=2 ) >= ntemp, &
      '2nd dimension of `selfen_dw` must be at least equal to length of array `temps`.' )
    call assert( size( selfen_dw, dim=1 ) >= nwfn, &
      '1st dimension of `selfen_dw` too small.' )

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
        if (ph_energy(imode, iq) <= eph_ph_energy_zero) cycle
        do is = 1, nspecies
          if (spmass(is) == 0.0_dp) cycle
          mi = 1.0_dp / (2.0_dp * ph_energy(imode, iq) * spmass(is))
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
          nq = bose_einstein( ph_energy(imode, iq), temps(itemp) )
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
            selfen_dw(ist, itemp) = selfen_dw(ist, itemp) - g(jst, ist, 1) / cmplx( de, eta, dp )
          end do
        end do
        !$omp end do
        !$omp end parallel
      end do
    end do
  end subroutine eph_else_gen_debye_waller_ahc
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
    call assert( size( selfen ) == n, &
      '`freqs` and `selfen` must be of equal size.' )
    call assert( size( sfun ) == n, &
      '`freqs` and `sfun` must be of equal size.' )

    ! infinitesimal to get peak even if Im(Sigma) is zero
    eps = sign( eps_zero, selfen(maxloc( abs( selfen%im ), dim=1 ))%im )

    sfun = piinv * (selfen%im + eps) / ((freqs - selfen%re)**2 + (selfen%im + eps)**2)
  end subroutine eph_else_gen_specfun
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! AUXILIARY PROCEDURES
  !
  !> Set default parameters for frequency grid used of electron self-energy calculation.
  pure subroutine eph_else_set_default_frequency_grid( fgrid )
    use eph_electrons, only: eph_el_energy_k

    use modinput, only: freq_grid_type
    !> frequency grid object
    type(freq_grid_type), intent(out) :: fgrid
  
    fgrid%type = 'density'          ! density based sampling
    fgrid%numpoints = 500           ! number of sampling points
    fgrid%range = [0.0_dp, 0.0_dp]  ! automatic range detection
    fgrid%padding = 0.1_dp          ! padding to add at both ends of range
    fgrid%lorentzwidth = 0.02_dp    ! width of Lorantzian density
  end subroutine eph_else_set_default_frequency_grid

  !> Get file name for binary file.
  function eph_else_get_binary_file_name( ngridbz, wfrange, nfreq, temps ) result( fname )
    !> Brillouin zone sampling used for integration
    integer, intent(in) :: ngridbz(3)
    !> Wannier function range
    integer, intent(in) :: wfrange(2)
    !> number of frequencies
    integer, intent(in) :: nfreq
    !> temperatures
    real(dp), intent(in) :: temps(:)
    !> file name
    character(:), allocatable :: fname

    integer :: ntemp, tmin, tmax, tstep
    character(8) :: bz_digit, wf_digit, freq_digit, temp_digit
    character(128) :: string

    ntemp = size(temps)
    tmin = nint( temps(1) )
    tmax = nint( temps(ntemp) )
    tstep = 0
    if (ntemp > 1) tstep = nint( abs( temps(ntemp) - temps(1) ) / (ntemp - 1) )

    ! set number of digits
    write( bz_digit, '(i8.8)' ) ceiling( log10( maxval( ngridbz ) + 1e-12_dp ) )
    write( wf_digit, '(i8.8)' ) ceiling( log10( maxval( wfrange ) + 1e-12_dp ) )
    write( freq_digit, '(i8.8)' ) ceiling( log10( nfreq + 1e-12_dp ) )
    write( temp_digit, '(i8.8)' ) ceiling( log10( maxval( [tmin, tmax, tstep] ) + 1e-12_dp ) )

    write( string, '("EPH_ELSE_",&
      &"BZ",3(i'//bz_digit//'.'//bz_digit//',"_"),&
      &"WF",2(i'//wf_digit//'.'//wf_digit//',"_"),&
      &"F",i'//freq_digit//'.'//freq_digit//',"_",&
      &"T",2(i'//temp_digit//'.'//temp_digit//',"_"),i'//temp_digit//'.'//temp_digit//',".OUT")' ) &
      ngridbz, wfrange, nfreq, tmin, tmax, tstep

    fname = trim( string )
  end function eph_else_get_binary_file_name

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
    real(dp), parameter :: sfun_wgt = 0.8_dp

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
        point_density = sfun_wgt * sfun_in + (1.0_dp - sfun_wgt) / (freqs_in(nfreq_in) - freqs_in(1))
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
        call terminate_if_false( .true., '(eph_else_resample) &
          Unsupported interpolation mode `'//method//'`.' )
    end select
    selfen_out = cmplx( fun1, fun2, dp )
  end subroutine eph_else_resample
  !-------------------------------------------------------------------------------- 

end module eph_electron_selfenergy
