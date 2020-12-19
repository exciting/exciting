!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! using the frequency convolution
!==================================================================
subroutine calcselfc_freqconv_cd(ikp,iq,mdim)
    use modinput
    use constants, only : zzero, pi
    use modmain, only : idxas, evalcr, efermi
    use modgw,   only : ibgw, nbgw, nstse, kset, kqset, freq, selfec, mwm, &
    &                   corind, freq_selfc, evalfv
    use mod_pade
    use mod_aaa_approximant
    ! input variables
    implicit none
    integer(4), intent(in) :: ikp
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: mdim
    ! local variables
    integer(4) :: ik, jk, jkp
    integer(4) :: ia, is, ias, ic, icg
    integer(4) :: ie1, ie2
    integer(4) :: iom, jom, jom1
    real(8)    :: enk, w, om1, om2
    complex(8) :: xnm(1:freq%nomeg), xnm_intp
    complex(8) :: sc, dfz, w_ac
    type(aaa_approximant) :: aaa

    ! k point
    ik = kset%ikp2ik(ikp)
    ! k-q point
    jk = kqset%kqid(ik,iq)
    jkp = kset%ik2ikp(jk)

    !------------------------
    ! loop over frequencies
    !------------------------

    do ie1 = ibgw, nbgw

      do ie2 = 1, mdim

        if ( ie2 <= nstse ) then
          !=============================
          ! Valence electron contribution
          !=============================
          enk = evalfv(ie2,jkp)-efermi
        else
          !=============================
          ! Core electron contribution
          !=============================
          icg = ie2-nstse
          is = corind(icg,1)
          ia = corind(icg,2)
          ic = corind(icg,3)
          ias = idxas(ia,is)
          enk = evalcr(ic,ias)-efermi
        end if ! val/cor

        xnm(:) = mwm(ie1,ie2,:)

        if (input%gw%selfenergy%actype == 'aaa') then
          call set_aaa_approximant(aaa, &
                                   cmplx(0.d0,freq%freqs,8), &
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

          om1 = 0.d0 ! Omega_{l}
          do jom = 1, freq%nomeg
            jom1 = min(jom+1, freq%nomeg) !  Omega_{l+1}
            om2 = 0.5d0 * (freq%freqs(jom1) + freq%freqs(jom))
            sc = sc + xnm(jom) * (atan(om2/(w-enk)) - atan(om1/(w-enk)))
            om1 = om2 ! next integration sub-interval
          end do
          sc = sc / pi

          !------------------------------------
          ! 2) contribution from W poles
          !------------------------------------
          w_ac = cmplx(abs(w-enk), 0.d0, 8) ! |w-e_nk|-i*eta

          ! Analytical continuation
          if (input%gw%selfenergy%actype == 'pade') then
            call pade_approximant(freq%nomeg, cmplx(0.d0,freq%freqs,8), xnm, &
                                  w_ac, xnm_intp, dfz)
          else if (input%gw%selfenergy%actype == 'aaa') then
            xnm_intp = get_aaa_approximant(aaa, w_ac)
            xnm_intp = conjg(xnm_intp)
          end if
          sc = sc + (theta(enk-w)*theta(-enk) - &
                     theta(w-enk)*theta(enk)) * xnm_intp

          ! sum over states
          selfec(ie1,iom,ikp) = selfec(ie1,iom,ikp) - sc

        end do ! frequency loop

        if (input%gw%selfenergy%actype == 'aaa') &
          call delete_aaa_approximant(aaa)

      end do ! ie2

    end do ! ie1

    return

contains

    real(8) function theta(x)
      real(8), intent(in) :: x
      if (x > 0.d0) then
        theta = 1.d0
      else
        theta = 0.d0
      end if
    end function

end subroutine
