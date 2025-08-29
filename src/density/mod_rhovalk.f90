module mod_rhovalk
  use asserts, only: assert
  use constants, only: real_zero, zone, zzero
  use mod_APW_LO, only: apwordmax, lofr, nlorb, lorbl
  use mod_atoms, only: natmtot, nspecies, natoms, idxas
  use mod_eigensystem, only: idxlo
  use mod_eigenvalue_occupancy, only: nstfv, nstsv
  use mod_gkvector, only: ngk, ngkmax, gkc, tpgkc, sfacgk
  use mod_Gvector, only: ngrtot
  use mod_kpoint, only: wkpt
  use mod_muffin_tin, only: nrcmtmax, lmmaxapw, idxlm, nrcmt, nrmt, lmmaxvr, nrmtmax
  use mod_SHT, only: zbshtvr, rfshtvr
  use mod_spin, only: ncmag, nspnfv, nspinor, ndmag
  use mod_timing, only: timerho
  use modinput, only: input, isspinspiral, issvlo
  use modmpi, only: mpi_env_k, distribute_loop
  use precision, only: dp, i32
  use svlo, only: get_num_of_basis_functions_sv

  implicit none

  private

  public :: rhovalk

  interface rhovalk
    module procedure :: rhovalk_spin_polarized
    module procedure :: rhovalk_non_spin_polarized
  end interface

contains

  !> Remaps non-spin-polarized wavefunctions to the spin polarized ones
  !> to be used as input arguments for `rhovalk_spin_polarized`
  subroutine rhovalk_non_spin_polarized ( ik, evecfv, occupations, rhomt, magmt, evecsv )
    !> k-point number
    integer(i32), intent (in) :: ik
    !> First-variational eigenvectors (nmatmax, nstfv)
    complex(dp), contiguous, target, intent (in) :: evecfv(:, :)
    !> State occupations (nstsv)
    real(dp), contiguous, intent(in) :: occupations(:)
    !> Muffin-tin charge density
    real(dp), contiguous, intent(inout) :: rhomt(:, :, :)
    !> Muffin-tin magnetisation vector field
    real(dp), contiguous, optional, intent(inout) :: magmt(:, :, :, :)
    !> Second-variational eigenvectors (nstfv, nstsv)
    complex(dp), contiguous, optional, intent(in) :: evecsv(:, :)

    integer(i32), parameter :: n_spin = 1
    complex(dp), contiguous, pointer :: ptr(:, :, :)

    ptr(1 : size( evecfv, 1 ), 1 : size( evecfv, 2 ), 1 : n_spin) => evecfv
    call rhovalk_spin_polarized( ik, ptr, occupations, rhomt, magmt, evecsv )

  end subroutine rhovalk_non_spin_polarized

  !> Generates the partial valence charge density from the eigenvectors at
  !> $k$-point {\tt ik}. In the muffin-tin region, the wavefunction is obtained
  !> in terms of its $(l,m)$-components from both the APW and local-orbital
  !> functions. Using a backward spherical harmonic transform (SHT), the
  !> wavefunction is converted to real-space and the density obtained from its
  !> modulus squared. This density is then transformed with a forward SHT and
  !> accumulated in the inout variable {\tt rhomt}.
  subroutine rhovalk_spin_polarized ( ik, evecfv, occupations, rhomt, magmt, evecsv )
    !> k-point number
    integer(i32), intent (in) :: ik
    !> First-variational eigenvectors (nmatmax, nstfv, n_spin)
    complex(dp), contiguous, intent (in) :: evecfv (:, :, :)
    !> State occupations (nstsv)
    real(dp), contiguous, intent(in) :: occupations(:)
    !> Muffin-tin charge density
    real(dp), contiguous, intent(inout) :: rhomt(:, :, :)
    !> Muffin-tin magnetisation vector field
    real(dp), contiguous, optional, intent(inout) :: magmt(:, :, :, :)
    !> Second-variational eigenvectors (nstfv, nstsv)
    complex(dp), contiguous, optional, intent(in) :: evecsv(:, :)

    integer(i32) :: nsd, ispn, jspn, is, ia, ias, ist
    integer(i32) :: ir, irc, itp, i, j, n, ilo, l, m, lm, nr
    real(dp) :: t1, ts0, ts1
    complex(dp) :: zt1, zt2, zt3
    integer(i32) :: num_of_basis_functions_sv, num_of_states
    logical, allocatable :: done(:, :)
    real(dp), allocatable :: rflm(:, :), rfmt(:, :, :)
    complex(dp), allocatable :: apwalm(:, :, :, :, :)
    complex(dp), allocatable :: wfmt1(:, :), wfmt2(:, :, :, :), wfmt3(:, :, :)
    real(dp) :: rhomt_k(lmmaxvr, nrmtmax, natmtot), magmt_k(lmmaxvr, nrmtmax, natmtot, ndmag)

    call timesec( ts0 )

    num_of_basis_functions_sv = get_num_of_basis_functions_sv()
    if ( input%groundstate%tevecsv ) then
      call assert( present( evecsv ), 'evecsv not present' )
      num_of_states = size( evecsv, 2 )
    else
      num_of_states = size( evecfv, 2 )
    end if

    if ( associated( input%groundstate%spin ) ) then
      call assert( present( magmt ), 'magmt not present' )
      magmt_k = real_zero
      if ( ncmag ) then
        nsd = 4
      else
        nsd = 2
      end if
    else
      nsd = 1
    end if

    rhomt_k = real_zero
    
    allocate( done(num_of_basis_functions_sv, nspnfv) )
    allocate( rflm(lmmaxvr, nsd) )
    allocate( rfmt(lmmaxvr, nrcmtmax, nsd) )
    allocate( apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv) )
    allocate( wfmt1(lmmaxvr, nrcmtmax) )
    if ( input%groundstate%tevecsv ) allocate( wfmt2(lmmaxvr, nrcmtmax, &
    & num_of_basis_functions_sv, nspnfv) )
    allocate( wfmt3(lmmaxvr, nrcmtmax, nspinor) )

    ! find the matching coefficients
    do ispn = 1, nspnfv
      call match( ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, &
        & ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn) )
    end do
    
    do is = 1, nspecies
      n = lmmaxvr * nrcmt(is)
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        done = .false.
        rfmt = real_zero
        do j = 1, num_of_states
          t1 = wkpt (ik) * occupations(j)
          if ( abs(t1) <= input%groundstate%epsocc ) cycle
          if ( input%groundstate%tevecsv ) then
            ! generate spinor wavefunction from second-variational eigenvectors
            wfmt3 = zzero
            do ispn = 1, nspinor
              if ( isspinspiral() ) then
                jspn = ispn
              else
                jspn = 1
              end if
              do ist = 1, nstfv
                i = (ispn - 1) * num_of_basis_functions_sv + ist
                zt1 = evecsv(i, j)
                if ( abs( real( zt1, dp ) ) + abs( aimag( zt1 ) ) > &
                  input%groundstate%epsocc ) then
                  if ( .not. done(ist, jspn) ) then
                    if ( issvlo() ) then
                      call wavefmt_apw( input%groundstate%lradstep, &
                        input%groundstate%lmaxvr, is, ia, ngk(jspn, ik), &
                        apwalm(:, :, :, :, jspn), evecfv(:, ist, jspn), lmmaxvr, wfmt1 )
                    else
                      call wavefmt( input%groundstate%lradstep, &
                        input%groundstate%lmaxvr, is, ia, ngk(jspn, ik), &
                        apwalm(:, :, :, :, jspn), evecfv(:, ist, jspn), lmmaxvr, wfmt1 )
                    end if
                    ! convert from spherical harmonics to spherical coordinates
                    call zgemm( 'N', 'N', lmmaxvr, nrcmt(is), lmmaxvr, zone, zbshtvr, &
                      lmmaxvr, wfmt1, lmmaxvr, zzero, wfmt2(:, :, ist, jspn), lmmaxvr )
                    done (ist, jspn) = .true.
                  end if
                  ! add to spinor wavefunction
                  call zaxpy( n, zt1, wfmt2(:, :, ist, jspn), 1, wfmt3(:, :, ispn), 1 )
                end if
              end do
              ! add local orbital contribution in case of a svlo calculation 
              if ( issvlo() ) then
                do ilo = 1, nlorb(is)
                  l = lorbl(ilo, is)
                  if ( l <= input%groundstate%lmaxvr ) then
                    do m= -l, l
                      lm = idxlm(l, m)
                      ist = nstfv + idxlo(lm, ilo, ias)
                      i = (ispn - 1) * num_of_basis_functions_sv + ist
                      zt1 = evecsv(i, j)
                      if ( abs( dble( zt1 ) ) + abs( aimag( zt1 ) ) > &
                        input%groundstate%epsocc ) then
                        if ( .not. done(ist, jspn) ) then
                          nr = 0
                          do ir = 1, nrmt(is), input%groundstate%lradstep
                            nr = nr + 1
                            wfmt2(:, nr, ist, jspn) = zbshtvr(:, lm) * lofr(ir, 1, ilo, ias)
                          end do
                          done(ist, jspn) = .true.
                        end if
                        call zaxpy( n, zt1, wfmt2(:, :, ist, jspn), 1, wfmt3(:, :, ispn), 1 )
                      end if
                    end do
                  end if
                end do
              end if
            end do
          else
            ! spin-unpolarised wavefunction
            call wavefmt ( input%groundstate%lradstep, &
              input%groundstate%lmaxvr, is, ia, ngk(1, ik), &
              apwalm, evecfv(:, j, 1), lmmaxvr, wfmt1 )
            ! convert from spherical harmonics to spherical coordinates
            call zgemm ( 'N', 'N', lmmaxvr, nrcmt(is), lmmaxvr, &
              zone, zbshtvr, lmmaxvr, wfmt1, lmmaxvr, zzero, wfmt3, lmmaxvr )
          end if
          
          ! add to the spin density matrix
          if ( associated( input%groundstate%spin ) ) then
            ! spin-polarised
            do irc = 1, nrcmt (is)
              do itp = 1, lmmaxvr
                zt1 = wfmt3(itp, irc, 1)
                zt2 = wfmt3(itp, irc, 2)
                zt3 = zt1 * conjg (zt2)
                rfmt(itp, irc, 1) = rfmt (itp, irc, 1) + t1 * (dble (zt1) ** 2 + aimag (zt1) ** 2)
                rfmt(itp, irc, 2) = rfmt (itp, irc, 2) + t1 * (dble (zt2) ** 2 + aimag (zt2) ** 2)
                if ( ncmag ) then
                  rfmt(itp, irc, 3) = rfmt(itp, irc, 3) + t1 * dble( zt3 )
                  rfmt (itp, irc, 4) = rfmt (itp, irc, 4) + t1 * aimag (zt3)
                end if
              end do
            end do
          else
            ! spin-unpolarised
            do irc = 1, nrcmt(is)
              do itp = 1, lmmaxvr
                zt1 = wfmt3(itp, irc, 1)
                rfmt(itp, irc, 1) = rfmt(itp, irc, 1) + t1 * (dble (zt1) ** 2 + aimag (zt1) ** 2)
              end do
            end do
          end if
        end do
        
        ! convert to spherical harmonics and add to rhomt_k and magmt_k
        irc = 0
        do ir = 1, nrmt(is), input%groundstate%lradstep
          irc = irc + 1
          do i = 1, nsd
            call dgemv( 'N', lmmaxvr, lmmaxvr, 1._dp, rfshtvr, &
              lmmaxvr, rfmt(:, irc, i), 1, 0._dp, rflm(:, i), 1 )
          end do
          if ( associated( input%groundstate%spin ) ) then
            ! spin-polarised
            if ( ncmag ) then
              magmt_k(:, ir, ias, 1) = magmt_k(:, ir, ias, 1) + 2._dp * rflm(:, 3)
              magmt_k(:, ir, ias, 2) = magmt_k(:, ir, ias, 2) - 2._dp * rflm(:, 4)
              magmt_k(:, ir, ias, 3) = magmt_k(:, ir, ias, 3) + rflm(:, 1) - rflm(:, 2)
            else
              magmt_k(:, ir, ias, 1) = magmt_k(:, ir, ias, 1) + rflm(:, 1) - rflm(:, 2)
            end if
            rhomt_k(:, ir, ias) = rhomt_k(:, ir, ias) + rflm(:, 1) + rflm(:, 2)
          else
            ! spin-unpolarised
            rhomt_k(:, ir, ias) = rhomt_k(:, ir, ias) + rflm(:, 1)
          end if
        end do
      end do
    end do

    rhomt = rhomt + rhomt_k
    if ( associated( input%groundstate%spin ) ) magmt = magmt + magmt_k
    
    call timesec (ts1)
    timerho = timerho + ts1 - ts0
  end subroutine rhovalk_spin_polarized

end module