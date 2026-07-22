!> Module to compute Eliashberg spectral function \(\alpha^2F_n({\bf k},\varepsilon,\omega)\) 
!> and related quantities.
module eph_eliashberg
  use eph_variables

  use precision, only: dp
#include "asserts.fpp"
  use modmpi

  implicit none
  private

  public :: eph_eliashberg_gen_a2F_tetrahedron, &
            eph_eliashberg_gen_phonon_coupling_tetrahedron

contains

  !================================================================================ 
  ! ELIASHBERG SPECTRAL FUNCTION
  !
  !> Calculate Eliashberg spectral function \(\alpha^2F^\pm_n({\bf k},\varepsilon,\omega)\) using
  !> tetrahedron integration.
  !>
  !> The Eliashberg spectral function is given by
  !> \[
  !>    \alpha^2F^\pm_n({\bf k},\varepsilon,\omega) = \sum\limits_{m,\nu} \int \frac{{\rm d}{\bf q}}{V_{\rm BZ}}
  !>    |g_{mn,\nu}({\bf k},{\bf q})|^2\, \delta(\omega - \omega_{\nu{\bf q}})\, \delta(\varepsilon \pm \omega_{\nu{\bf q}} - \epsilon_{m{\bf k}+{\bf q}}) \;,
  !> \]
  !> where the BZ integral is computed using tetrahedron integration.
  !> The sign \(\pm\) corresponds to phonon absorption/emission. If `sig = 0`, the quasielastic approximation
  !> is applied, i.e., the phonon frequencies \(\pm \omega_{\nu{\bf q}}\) in the electron Dirac delta 
  !> are neglected.
  !>
  !> This subroutine does not include any internal MPI parallelism. But the sums over states \(m\) and phonon modes \(\nu\) can be 
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_gen_electron_coupling_strength(subroutine)]].
  subroutine eph_eliashberg_gen_a2F_tetrahedron( el_freqs, ph_freqs, el_energy_k, el_energy_kq, ph_energy_q, ephmat, tset, a2F, &
      sig )
    use mod_opt_tetra, only: t_set, opt_tetra_int_dbldelta
    use math_utils, only: get_degeneracies
    !> electron frequencies \(\varepsilon\)
    real(dp), intent(in) :: el_freqs(:)
    !> phonon frequencies \(\omega\)
    real(dp), intent(in) :: ph_freqs(:)
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
    !> Eliashberg spectral function \(\alpha^2F^\pm_n({\bf k},\varepsilon,\omega)\)
    real(dp), intent(out) :: a2F(:, :, :)
    !> `> 0` for phonon absorption   
    !> `< 0` for phonon emission   
    !> `= 0` for quasielastic approximation (default)
    integer, optional, intent(in) :: sig

    integer :: nwfk, nwfkq, nmode, nfreqel, nfreqph, nkq, iq, imode, ist, jst, signum, i
 
    integer, allocatable :: deg(:,:)
    real(dp), allocatable :: e(:,:), matr(:,:,:,:), intr(:,:,:,:,:)

    ! set defaults
    signum = 0
    if (present(sig)) signum = sig

    ! set matrix sizes
    nwfk = size( ephmat, dim=2 )
    nwfkq = size( ephmat, dim=1 )
    nmode = size( ephmat, dim=3 )
    nkq = size( ephmat, dim=4 )
    nfreqel = size( el_freqs )
    nfreqph = size( ph_freqs )

    ! check input
    CALL_ASSERT( size( el_energy_k ) == nwfk, '1st dimension of `el_energy_k` must equal 2nd dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=1 ) == nwfkq, '1st dimension of `el_energy_kq` must equal 1st dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=2 ) == nkq, '2nd dimension of `el_energy_kq` must equal 4th dimension of `ephmat`.' )
    CALL_ASSERT( size( ph_energy_q, dim=1 ) == nmode, '1st dimension of `ph_energy_q` must equal 3rd dimension of `ephmat`.' )
    CALL_ASSERT( size( ph_energy_q, dim=2 ) == nkq, '2nd dimension of `ph_energy_q` must equal number of points `qset`.' )
    CALL_ASSERT( size( a2F, dim=1 ) >= nfreqph, '1st dimension of `a2F` must be at least equal to length of array `ph_freqs`.' )
    CALL_ASSERT( size( a2F, dim=2 ) >= nfreqel, '2nd dimension of `a2F` must be at least equal to length of array `el_freqs`.' )
    CALL_ASSERT( size( a2F, dim=3 ) >= nwfk, '3rd dimension of `a2F` too small.' )

    a2F(:nfreqph, :nfreqel, :nwfk) = 0.0_dp

    allocate( e(nwfkq, nkq), matr(nwfk, nwfkq, 1, nkq) )
    allocate( intr(nfreqel, nfreqph, nwfk, nwfkq, 1) )

    deg = get_degeneracies( el_energy_k, eph_el_degtol )

    do imode = 1, nmode
      !$omp parallel default( shared )
      !$omp do
      do iq = 1, nkq
        if (signum > 0) then
          ! Note: delta(epsilon - (epsilon_{m k+q} - omega_{nu q})), hence '-'
          e(:, iq) = el_energy_kq(:, iq) - ph_energy_q(imode, iq)
        else if (signum < 0) then
          ! Note: delta(epsilon - (epsilon_{m k+q} + omega_{nu q})), hence '+'
          e(:, iq) = el_energy_kq(:, iq) + ph_energy_q(imode, iq)
        else
          e(:, iq) = el_energy_kq(:, iq)
        end if
        do jst = 1, nwfkq
          matr(:, jst, 1, iq) = ephmat(jst, :, imode, iq)%re**2 + ephmat(jst, :, imode, iq)%im**2
          ! average over degenerate states
          do i = 1, size(deg, dim=2)
            matr(deg(1, i):deg(2, i), jst, 1, iq) = sum( matr(deg(1, i):deg(2, i), jst, 1, iq) ) / deg(3, i)
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
      call opt_tetra_int_dbldelta( tset, nkq, nwfkq, e, 1, ph_energy_q(imode:imode, :), nfreqel, el_freqs, nfreqph, ph_freqs, &
        nwfk, intr, fun=matr )
      do jst = 1, nwfkq
        do ist = 1, nwfk
          a2F(:nfreqph, :nfreqel, ist) = a2F(:nfreqph, :nfreqel, ist) + transpose( intr(:, :, ist, jst, 1) )
        end do
      end do
    end do

    deallocate( e, matr, intr )
  end subroutine eph_eliashberg_gen_a2F_tetrahedron

  !================================================================================ 
  ! PHONON EPH COUPLING STRENGTH
  !
  !> Calculate phonon resolved EPH coupling strength \(\lambda^\pm_{\nu{\bf q}}\) using
  !> tetrahedron integration.
  !>
  !> The phonon EPH coupling strength is given by
  !> \[
  !>    \lambda^\pm_{\nu{\bf q}}(\varepsilon) = \frac{1}{N(\varepsilon)\, \omega_{\nu{\bf q}}} \sum\limits_{m,n} \int \frac{{\rm d}{\bf k}}{V_{\rm BZ}} 
  !>    |g_{mn,\nu}({\bf k},{\bf q})|^2\, \delta(\varepsilon - \epsilon_{n{\bf k}})\, \delta(\epsilon_{n{\bf k}} \pm \omega_{\nu{\bf q}} - \epsilon_{m{\bf k}+{\bf q}}) \;,
  !> \]
  !> where the BZ integral is computed using tetrahedron integration.
  !> The sign \(\pm\) corresponds to phonon absorption/emission. If `sig = 0`, the quasielastic approximation
  !> is applied, i.e., the phonon frequencies \(\pm \omega_{\nu{\bf q}}\) in the electron Dirac delta 
  !> are neglected.
  !>
  !> This subroutine does not include any internal MPI parallelism. But the sums over states \(m\) and phonon modes \(\nu\) can be 
  !> be distributed by only passing corresponding subsets. See e.g. [[eph_gen_electron_coupling_strength(subroutine)]].
  subroutine eph_eliashberg_gen_phonon_coupling_tetrahedron( el_freqs, el_energy_k, el_energy_kq, ph_energy_q, ephmat, tset, lambda, &
      sig )
    use mod_opt_tetra, only: t_set, opt_tetra_int_delta, opt_tetra_int_dbldelta
    use math_utils, only: get_degeneracies
    !> electron frequencies \(\varepsilon\)
    real(dp), intent(in) :: el_freqs(:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: el_energy_k(:, :)
    !> electron energies \(\epsilon_{m{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: el_energy_kq(:, :)
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy_q(:)
    !> EPH matrix elements \(g_{mn,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: ephmat(:, :, :, :)
    !> set of tetrahedra for integration
    type(t_set), intent(in) :: tset
    !> Eliashberg spectral function \(\lambda^\pm_{\nu{\bf q}}(\varepsilon)\)
    real(dp), intent(out) :: lambda(:, :)
    !> `> 0` for phonon absorption   
    !> `< 0` for phonon emission   
    !> `= 0` for quasielastic approximation (default)
    integer, optional, intent(in) :: sig

    integer :: nwfk, nwfkq, nmode, nfreqel, nk, ik, imode, ist, jst, signum, i
 
    integer, allocatable :: deg(:,:)
    real(dp), allocatable :: e(:,:), dos(:,:,:), matr(:,:,:,:), intr(:,:,:,:,:)

    ! set defaults
    signum = 0
    if (present(sig)) signum = sig

    ! set matrix sizes
    nwfk = size( ephmat, dim=2 )
    nwfkq = size( ephmat, dim=1 )
    nmode = size( ephmat, dim=3 )
    nk = size( ephmat, dim=4 )
    nfreqel = size( el_freqs )

    ! check input
    CALL_ASSERT( size( el_energy_k, dim=1 ) == nwfk, '1st dimension of `el_energy_k` must equal 2nd dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=1 ) == nwfkq, '1st dimension of `el_energy_kq` must equal 1st dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_k, dim=2 ) == nk, '2nd dimension of `el_energy_k` must equal 4th dimension of `ephmat`.' )
    CALL_ASSERT( size( el_energy_kq, dim=2 ) == nk, '2nd dimension of `el_energy_kq` must equal 4th dimension of `ephmat`.' )
    CALL_ASSERT( size( ph_energy_q ) == nmode, '1st dimension of `ph_energy_q` must equal 3rd dimension of `ephmat`.' )
    CALL_ASSERT( size( lambda, dim=1 ) >= nfreqel, '1st dimension of `lambda` must be at least equal to length of array `el_freqs`.' )
    CALL_ASSERT( size( lambda, dim=2 ) >= nmode, '2nd dimension of `lambda` too small.' )

    lambda(:nfreqel, :nmode) = 0.0_dp

    allocate( e(nwfkq, nk), matr(nmode, nwfkq, 1, nk) )
    if (signum == 0) then
      allocate( intr(1, nfreqel, nmode, nwfkq, 1) )
    else
      allocate( intr(nmode, nfreqel, nmode, nwfkq, 1) )
    end if

    ! compute DOS
    allocate( dos(nfreqel, 1, nwfk) )
    call opt_tetra_int_delta( tset, nk, nwfk, el_energy_k, nfreqel, el_freqs, 1, dos )
    dos(:, :, 1) = sum( dos, dim=3 )

    deg = get_degeneracies( ph_energy_q, eph_ph_degtol )
    do ist = 1, nwfk
      !$omp parallel default( shared )
      !$omp do
      do ik = 1, nk
        e(:, ik) = el_energy_kq(:, ik) - el_energy_k(ist, ik)
        do jst = 1, nwfkq
          matr(:, jst, 1, ik) = ephmat(jst, ist, :, ik)%re**2 + ephmat(jst, ist, :, ik)%im**2
          ! average over degenerate phonon modes
          do i = 1, size(deg, dim=2)
            matr(deg(1, i):deg(2, i), jst, 1, ik) = sum( matr(deg(1, i):deg(2, i), jst, 1, ik) ) / deg(3, i)
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
      if (signum > 0) then
        ! Note: delta(epsilon_{m k+q} - epsilon_{n k} - omega_{nu q}), hence '-'
        call opt_tetra_int_dbldelta( tset, nk, nwfkq, e, 1, el_energy_k(ist:ist, :), nmode, -ph_energy_q, nfreqel, el_freqs, &
          nmode, intr, fun=matr )
      else if (signum < 0) then
        ! Note: delta(epsilon_{m k+q} - epsilon_{n k} + omega_{nu q}), hence '+'
        call opt_tetra_int_dbldelta( tset, nk, nwfkq, e, 1, el_energy_k(ist:ist, :), nmode,  ph_energy_q, nfreqel, el_freqs, &
          nmode, intr, fun=matr )
      else
        call opt_tetra_int_dbldelta( tset, nk, nwfkq, e, 1, el_energy_k(ist:ist, :), 1, [0.0_dp], nfreqel, el_freqs, &
          nmode, intr, fun=matr )
      end if
      do imode = 1, nmode
        do jst = 1, nwfkq
          if (signum == 0) then
            lambda(:nfreqel, imode) = lambda(:nfreqel, imode) + intr(1, :, imode, jst, 1)
          else
            lambda(:nfreqel, imode) = lambda(:nfreqel, imode) + intr(imode, :, imode, jst, 1)
          end if
        end do
      end do
    end do
    do imode = 1, nmode
      lambda(:nfreqel, imode) = dble( lambda(:nfreqel, imode) / cmplx( dos(:, 1, 1) * ph_energy_q(imode), eph_ph_degtol, dp) ) 
    end do

    deallocate( e, matr, intr )
  end subroutine eph_eliashberg_gen_phonon_coupling_tetrahedron
  !-------------------------------------------------------------------------------- 

end module eph_eliashberg
