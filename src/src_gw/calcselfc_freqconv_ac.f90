!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! using the frequency convolution
!==================================================================
subroutine calcselfc_freqconv_ac(ikp,iq,mdim)
    use modinput
    use modmain, only : idxas, evalcr, efermi
    use constants, only : zzero, pi
    use modgw,   only : ibgw, nbgw, nstse, kset, kqset, freq, selfec, mwm, &
    &                   ncg, corind, fdebug, freq_selfc, evalfv
    ! input variables
    implicit none
    integer(4), intent(in) :: ikp
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: mdim
    ! local variables
    integer(4) :: ik, jk, jkp
    integer(4) :: ia, is, ias, ic, icg
    integer(4) :: ie1, ie2, i1, i2, n, m
    integer(4) :: iom, jom, kom
    real(8)    :: enk, wdiff, w_sc, f1, f2
    complex(8) :: xnm(1:freq%nomeg)
    complex(8) :: sc, zt1, zt2

    ! k point
    ik = kset%ikp2ik(ikp)
    ! k-q point
    jk = kqset%kqid(ik,iq)
    jkp = kset%ik2ikp(jk)

    !------------------------
    ! loop over frequencies
    !------------------------

    do ie1 = ibgw, nbgw

      ! sum over states
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

        ! Re W
        xnm(:) = mwm(ie1,ie2,:)
        ! xnm(:) = cmplx( dble(mwm(ie1,ie2,:)), 0.d0, 8)

        ! for each frequency
        do iom = 1, freq_selfc%nomeg

          w_sc = freq_selfc%freqs(iom)

          !--------------------------------
          ! frequency convolution integral
          !--------------------------------

          ! (enk-iu)
          zt1 = cmplx( enk, -w_sc, 8)

          sc = zzero
          do jom = 1, freq%nomeg
            zt2 = freq%womeg(jom) / ( freq%freqs(jom)**2 + zt1**2 )
            sc = sc + (xnm(jom)-xnm(iom)) * zt2
          end do
          sc = sc*zt1/pi + xnm(iom)*sign(0.5d0,enk)

          ! sum over states
          selfec(ie1,iom,ikp) = selfec(ie1,iom,ikp) + sc

        end do ! frequency loop

      end do ! ie2

    end do ! ie1

    return
end subroutine
