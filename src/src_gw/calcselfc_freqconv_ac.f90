!> Calculate the q-dependent correlation self-energy term using frequency convolution.
subroutine calcselfc_freqconv_ac(ikp, iq)
    use modinput, only: input
    use mod_atoms, only: idxas
    use evgw0_validation, only: has_checkpoint_qp_energy_for_band
    use mod_eigenvalue_occupancy, only: efermi
    use mod_corestate, only: evalcr
    use constants, only : zzero, pi
    use modgw,   only : ibgw, nbgw, kset, kqset, freq, fdebug
    use mod_selfenergy, only: selfec, mwm, freq_selfc
    use mod_core_states, only: ncg, corind
    use mod_bands, only: nstse, evalfv
#include "mod_gw_degeneracies.inc"
    use mod_gw_degeneracies, only: get_degenerate_limits_qp_interval_ikp, &
                                   degenerate_subspaces
    use precision, only: i32, dp
    use self_consistent_eigenvalue_gw0, only: use_evgw0_input_qp, get_evalqp_evgw0_pointer

    implicit none

    !> Irreducible k-point index.
    integer(i32), intent(in) :: ikp
    !> q-point index.
    integer(i32), intent(in) :: iq

    integer(i32) :: ik, jk, jkp
    integer(i32) :: ia, is, ias, ic, icg
    integer(i32) :: ie1, ie2
    integer(i32) :: iom, jom
    real(dp)     :: enk, wdiff, w_sc
    complex(dp)  :: xnm(1:freq%nomeg)
    complex(dp)  :: sc, zt1, zt2
    ! For the averaging over degenerate states
    integer(i32) :: ispace_init, ispace_final, ispace, lowband, upband, size_deg
    real(dp), pointer :: evalqp_evgw0(:, :)

    if ( use_evgw0_input_qp() ) then
      call get_evalqp_evgw0_pointer(evalqp_evgw0)
    end if

    ! k point
    ik = kset%ikp2ik(ikp)
    ! k-q point
    jk = kqset%kqid(ik,iq)
    jkp = kset%ik2ikp(jk)

    ! First obtain the limits for degenerate subspaces for the given irreducible point. 
    call get_degenerate_limits_qp_interval_ikp(ikp, ispace_init, ispace_final)

    !-------------------------------
    ! Loop over degenerate subspaces
    !-------------------------------
    do ispace = ispace_init, ispace_final

      lowband  = degenerate_subspaces(1, ispace, ikp)
      upband   = degenerate_subspaces(2, ispace, ikp)
      size_deg = degenerate_subspaces(3, ispace, ikp)

      ! Sum over states in the degenerate subspace
      do ie1 = lowband, upband
        
        do ie2 = lbound( mwm, 2 ), ubound( mwm, 2 )

          if ( ie2 <= nstse ) then
            !=============================
            ! Valence electron contribution
            !=============================
            if ( use_evgw0_input_qp() ) then
              if ( has_checkpoint_qp_energy_for_band(ie2, lbound(evalqp_evgw0, 1), ubound(evalqp_evgw0, 1)) ) then
                enk = evalqp_evgw0(ie2,jkp)
              else
                enk = evalfv(ie2,jkp)
              end if
            else
              enk = evalfv(ie2,jkp)
            end if
          else
            !=============================
            ! Core electron contribution
            !=============================
            icg = ie2-nstse
            is = corind(icg,1)
            ia = corind(icg,2)
            ic = corind(icg,6)
            ias = idxas(ia,is)
            enk = evalcr(ic,ias) - efermi
          end if ! val/cor

          xnm(:) = mwm(ie1,ie2,:)

          ! for each frequency
          do iom = lbound( selfec, 2 ), ubound( selfec, 2 )

            w_sc = freq_selfc%freqs(iom)

            !--------------------------------
            ! frequency convolution integral
            !--------------------------------

            ! (enk-iu)
            zt1 = cmplx( enk, -w_sc, dp)

            sc = zzero
            do jom = 1, freq%nomeg
              zt2 = freq%womeg(jom) / ( freq%freqs(jom)**2 + zt1**2 )
              sc = sc + (xnm(jom)-xnm(iom)) * zt2
            end do
            sc = sc*zt1/pi + xnm(iom) * sign(0.5_dp,enk)

            ! Add the contribution to the degenerate subspace
            ! Note that the range macro "QP_ADJUST_RANGE" is found in mod_gw_degeneracies.inc
            selfec(QP_ADJUST_RANGE(lowband,upband),iom,ikp) = & 
                selfec(QP_ADJUST_RANGE(lowband,upband),iom,ikp) + sc / size_deg

          end do ! frequency loop

        end do ! ie2
      end do ! ie1

    end do ! ispace

    return
end subroutine
