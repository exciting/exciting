
subroutine calcselfc_freqconv_v2(ikp, iq, mdim)
    use modinput
    use modmain, only : pi, zzero, zi, evalsv, idxas, evalcr, efermi
    use modgw,   only : ibgw, nbgw, nstse, kset, kqset, freq, selfec, mwm, &
                        ncg, corind, fdebug, freq_selfc
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
    real(8)    :: enk, om1, om2
    complex(8) :: sc, w, zt1, zt2

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
          enk = evalsv(ie2,jkp)-efermi
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
          
        ! Self-energy frequency grid
        do iom = 1, freq_selfc%nomeg

          !------------------------
          ! Frequency convolution
          !------------------------

          ! (enk-iu)
          w = enk - zi*freq_selfc%freqs(iom)

          ! Omega_{l}
          om1 = 0.d0 ! freq%freqs(1)
          
          sc = zzero
          do jom = 1, freq%nomeg ! convolution integral

            !  Omega_{l+1}
            jom1 = min(jom+1, freq%nomeg)
            om2 = 0.5d0 * (freq%freqs(jom1) + freq%freqs(jom))
            ! print*, 'jom=', jom, freq%freqs(jom), om1, om2

            zt1 = om1 / w
            zt2 = om2 / w
            sc = sc + mwm(ie1,ie2,jom) * (arctan(zt2) - arctan(zt1))

            ! next integration point: Omega_{l}
            om1 = om2

          end do

          ! sum over states
          selfec(ie1,iom,ikp) = selfec(ie1,iom,ikp) + sc / pi

        end do ! frequency loop

      end do ! ie2

    end do ! ie1

    if (input%gw%debug) then
      write(fdebug,*) 'CORRELATION SELF-ENERGY: iq=', iq, ' ikp=', ikp
      write(fdebug,*) 'omega    state    Sigma_c'
      do iom = 1, freq%nomeg
        do ie1 = ibgw, nbgw
          write(fdebug,*) iom, ie1, selfec(ie1,iom,ikp)
        end do
      end do
      write(fdebug,*)
    end if
    
    return

contains

    complex(8) function arctan(z)
      implicit none
      complex(8) :: z
      complex(8), parameter :: zi = cmplx(0.d0, 1.d0, 8)
      arctan = 0.5d0 * zi * ( log(1.d0-zi*z) - log(1.d0+zi*z) )
    end function

end subroutine
