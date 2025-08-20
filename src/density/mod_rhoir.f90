module mod_rhoir

  use asserts, only: assert
  use constants, only: real_zero, zzero
  use m_zfftifc, only: zfftifc
  use mod_eigenvalue_occupancy, only: nstfv, nstsv
  use mod_gkvector, only: ngk, igkig
  use mod_Gvector, only: ngrtot, igfft, ngrid
  use mod_kpoint, only: wkpt
  use mod_lattice, only: omega
  use mod_spin, only: ncmag, nspinor, ndmag
  use mod_timing, only: timerho
  use modinput, only: input, isspinspiral
  use modmpi, only: mpi_env_k, distribute_loop
  use precision, only: dp, i32
  use svlo, only: get_num_of_basis_functions_sv

  implicit none

  private

  interface genrhoir
    module procedure :: genrhoir_spin_polarized
    module procedure :: genrhoir_non_spin_polarized
  end interface

  public :: genrhoir

contains

  !> Remaps non-spin-polarized wavefunctions to the spin polarized ones
  !> to be used as input arguments for `genrhoir_spin_polarized`
  subroutine genrhoir_non_spin_polarized ( ik, evecfv, occupations, rhoir, magir, evecsv )
    !> k-point number
    integer(i32), intent (in) :: ik
    !> First-variational eigenvectors (nmatmax, nstfv)
    complex(dp), contiguous, target, intent (in) :: evecfv(:, :)
    !> State occupations (nstsv)
    real(dp), contiguous, intent(in) :: occupations(:)
    !> Interstitial charge density
    real(dp), contiguous, intent(inout) :: rhoir(:)
    !> Interstitial magnetisation vector field
    real(dp), contiguous, optional, intent(inout) :: magir(:, :)
    !> Second-variational eigenvectors (nstfv, nstsv)
    complex(dp), contiguous, optional, intent(in) :: evecsv(:, :)

    integer(i32), parameter :: n_spin = 1
    complex(dp), contiguous, pointer :: ptr(:, :, :)

    ptr(1 : size( evecfv, 1 ), 1 : size( evecfv, 2 ), 1 : n_spin) => evecfv
    call genrhoir_spin_polarized( ik, ptr, occupations, rhoir, magir, evecsv )

  end subroutine genrhoir_non_spin_polarized

  !> Generates the partial valence charge density from the eigenvectors at
  !> k$-point {\tt ik}. The wavefunction in real-space is
  !> obtained from a Fourier transform of the sum of APW functions. The
  !> interstitial density is added to the inout array {\tt rhoir}. See routines
  !> {\tt wavefmt}, {\tt genshtmat} and {\tt seceqn}.
  subroutine genrhoir_spin_polarized ( ik, evecfv, occupations, rhoir, magir, evecsv )
    !> k-point number
    integer(i32), intent(in) :: ik
    !> First-variational eigenvectors (nmatmax, nstfv, n_spin)
    complex(dp), contiguous, intent (in) :: evecfv (:, :, :)
    !> State occupations (nstsv)
    real(dp), contiguous, intent(in) :: occupations(:)
    !> Interstitial charge density
    real(dp), contiguous, intent(inout) :: rhoir(:)
    !> Interstitial magnetisation vector field
    real(dp), contiguous, optional, intent(inout) :: magir(:, :)
    !> Second-variational eigenvectors (nstfv, nstsv)
    complex(dp), contiguous, optional, intent(in) :: evecsv(:, :)

    integer(i32) :: nsd, ispn, jspn, ist, ir, igk, ifg, i, j
    real(dp) :: t1, t2, t3, t4, ts0, ts1
    complex(dp) :: zt1, zt2, zt3
    integer(i32) :: num_of_basis_functions_sv, num_of_states
    complex(dp), allocatable :: zfft(:, :)
    real(dp) :: rhoir_k(ngrtot), magir_k(ngrtot, ndmag)

    call timesec( ts0 )
  
    num_of_basis_functions_sv = get_num_of_basis_functions_sv()
    if ( input%groundstate%tevecsv ) then
      call assert( present( evecsv ), 'evecsv not present' )
      num_of_states = size( evecsv, 2 )
    else
      num_of_states = size( evecfv, 2 )
    end if

    if ( associated( input%groundstate%spin ) ) then
      call assert( present( magir ), 'magir not present' )
      magir_k = real_zero
      if ( ncmag ) then
        nsd = 4
      else
        nsd = 2
      end if
    else
      nsd = 1
    end if

    rhoir_k = real_zero
    allocate ( zfft(ngrtot, nspinor) )
    
    do j = 1, num_of_states
      t1 = wkpt (ik) * occupations(j)
      if ( abs( t1 ) <= input%groundstate%epsocc ) cycle
      
      t2 = t1 / omega
      zfft = zzero
      if ( input%groundstate%tevecsv ) then
        ! generate spinor wavefunction from second-variational eigenvectors
        do ispn = 1, nspinor
          if ( isspinspiral() ) then
            jspn = ispn
          else
            jspn = 1
          end if
          do ist = 1, nstfv
            i = (ispn - 1) * num_of_basis_functions_sv + ist
            zt1 = evecsv (i, j)
            if ( abs( real( zt1, dp ) ) + abs( aimag( zt1 ) ) > &
              input%groundstate%epsocc ) then
              do igk = 1, ngk(jspn, ik)
                ifg = igfft(igkig(igk, jspn, ik))
                zfft(ifg, ispn) = zfft(ifg, ispn) + zt1 * evecfv(igk, ist, jspn)
              end do
            end if
          end do
        end do
      else
        ! spin-unpolarised wavefunction
        do igk = 1, ngk(1, ik)
          ifg = igfft(igkig(igk, 1, ik))
          zfft(ifg, 1) = evecfv(igk, j, 1)
        end do
      end if

      ! Fourier transform wavefunction to real-space
      do ispn = 1, nspinor
        call zfftifc( 3, ngrid, 1, zfft(:, ispn) )
      end do
      
      if ( associated( input%groundstate%spin ) ) then
        ! spin-polarised
        do ir = 1, ngrtot
          zt1 = zfft (ir, 1)
          zt2 = zfft (ir, 2)
          zt3 = zt1 * conjg(zt2)
          t3 = dble (zt1) ** 2 + aimag (zt1) ** 2
          t4 = dble (zt2) ** 2 + aimag (zt2) ** 2
          rhoir_k(ir) = rhoir_k(ir) + t2 * (t3 + t4)
          if ( ncmag ) then
            magir_k(ir, 1) = magir_k(ir, 1) + 2._dp * t2 * dble( zt3 )
            magir_k(ir, 2) = magir_k(ir, 2) - 2._dp * t2 * aimag( zt3 )
            magir_k(ir, 3) = magir_k(ir, 3) + t2 * (t3 - t4)
          else
            magir_k(ir, 1) = magir_k(ir, 1) + t2 * (t3 - t4)
          end if
        end do
      else
        ! spin-unpolarised
        do ir = 1, ngrtot
          zt1 = zfft(ir, 1)
          rhoir_k(ir) = rhoir_k(ir) + t2 * ( dble (zt1) ** 2 + aimag (zt1) ** 2 )
        end do
      end if
    end do
  
    rhoir = rhoir + rhoir_k
    if ( associated( input%groundstate%spin ) ) magir = magir + magir_k

    call timesec( ts1 )
    timerho = timerho + ts1 - ts0
  end subroutine

end module
