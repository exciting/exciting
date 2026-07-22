!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! using the frequency convolution
!==================================================================
subroutine calcselfc_freqconv_cd(ikp,iq,mdim)
    use modinput, only: input
    use constants, only : zzero, pi, real_one, real_zero
    use mod_atoms, only: idxas
    use mod_eigenvalue_occupancy, only: efermi
    use mod_corestate, only: evalcr
    use modgw,   only : ibgw, nbgw, kset, kqset, freq, fdebug
    use mod_selfenergy, only: selfec, mwm, freq_selfc
    use mod_core_states, only: ncg, corind
    use mod_bands, only: nstse, evalfv
#include "mod_gw_degeneracies.inc"
    use mod_gw_degeneracies, only: get_degenerate_limits_qp_interval_ikp, &
                                   degenerate_subspaces
    use mod_aaa_approximant, only: aaa_approximant, set_aaa_approximant, &
                                   init_aaa_approximant, delete_aaa_approximant, &
                                   get_aaa_approximant
    use mod_pade, only: pade_approximant
    use precision, only: i32, dp
    ! input variables
    implicit none
    integer(i32), intent(in) :: ikp
    integer(i32), intent(in) :: iq
    integer(i32), intent(in) :: mdim
    ! local variables
    integer(i32) :: ia, is, ias, ic, icg
    integer(i32) :: ik, jk, jkp
    integer(i32) :: ie1, ie2
    integer(i32) :: iom, jom, jom1
    real(dp)    :: enk, w, om1, om2
    complex(dp) :: xnm(1:freq%nomeg), xnm_intp
    complex(dp) :: sc, dfz, w_ac
    type(aaa_approximant) :: aaa
    ! For the averaging over degenerate states
    integer(i32) :: ispace_init, ispace_final, ispace, lowband, upband, size_deg

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

        do ie2 = 1, mdim

          if ( ie2 <= nstse ) then
            !=============================
            ! Valence electron contribution
            !=============================
            enk = evalfv(ie2,jkp)
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

          if (input%gw%selfenergy%actype == 'aaa') then
            call set_aaa_approximant(aaa, &
                                    cmplx(real_zero,freq%freqs,dp), &
                                    xnm, &
                                    input%gw%selfenergy%tol)
          end if

          ! Self-energy frequency grid
          do iom = 1, freq_selfc%nomeg

            w = freq_selfc%freqs(iom)

            !------------------------------------
            ! 1) frequency convolution integral
            !------------------------------------
            sc = zzero

            om1 = real_zero ! Omega_{l}
            do jom = 1, freq%nomeg
              jom1 = min(jom+1, freq%nomeg) !  Omega_{l+1}
              om2 = 0.5_dp * (freq%freqs(jom1) + freq%freqs(jom))
              sc = sc + xnm(jom) * (atan(om2/(w-enk)) - atan(om1/(w-enk)))
              om1 = om2 ! next integration sub-interval
            end do
            sc = sc / pi

            !------------------------------------
            ! 2) contribution from W poles
            !------------------------------------
            w_ac = cmplx(abs(w-enk), real_zero, dp) ! |w-e_nk|-i*eta

            ! Analytical continuation
            if (input%gw%selfenergy%actype == 'pade') then
              call pade_approximant(freq%nomeg, cmplx(real_zero,freq%freqs,dp), xnm, &
                                    w_ac, xnm_intp, dfz)
            else if (input%gw%selfenergy%actype == 'aaa') then
              xnm_intp = get_aaa_approximant(aaa, w_ac)
              xnm_intp = conjg(xnm_intp)
            end if
            sc = sc + (theta(enk-w)*theta(-enk) - &
                      theta(w-enk)*theta(enk)) * xnm_intp

            ! sum over states
            selfec(QP_ADJUST_RANGE(lowband,upband),iom,ikp) =  &
              selfec(QP_ADJUST_RANGE(lowband,upband),iom,ikp) - sc / size_deg

          end do ! frequency loop

          if (input%gw%selfenergy%actype == 'aaa') &
            call delete_aaa_approximant(aaa)

        end do ! ie2

      end do ! ie1

    end do ! ispace

    return

contains

    real(dp) function theta(x)
      real(dp), intent(in) :: x
      if (x > real_zero) then
        theta = real_one
      else
        theta = real_zero
      end if
    end function

end subroutine
