!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! using the frequency convolution
!==================================================================
subroutine calcselfc_freqconv(ikp,iq,mdim)
    use modinput
    use modmain, only : pi, zzero, evalsv, idxas, evalcr, efermi
    use modgw,   only : ibgw, nbgw, nstse, kset, kqset, freq, selfec, mwm, &
    &                   ncg, corind, fdebug, freq_selfc
    use mod_aaa_approximant
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
    complex(8) :: xnm(1:freq%nomeg), xnm_intp
    complex(8) :: sc, zt1, zt2

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
          
        xnm(:) = mwm(ie1,ie2,:)

        if (ikp==2 .and. iq==2 .and. ie1==1 .and. ie2==1) then
          do iom = 1, freq%nomeg
            write(10,'(3f18.6)') freq%freqs(iom), xnm(iom)
          end do
        end if

        ! AAA approximant
        call set_aaa_approximant(aaa, cmplx(freq%freqs, 0.d0, 8), xnm)

        ! Self-energy frequency grid
        do iom = 1, freq_selfc%nomeg

          w_sc = freq_selfc%freqs(iom)

          !----------------------
          ! Interpolation block
          !----------------------

          ! Linear interpolation
          ! if ( w_sc < freq%freqs(1) ) then
          !   kom = 1
          ! else if ( w_sc > freq%freqs(freq%nomeg) )  then
          !   kom = freq%nomeg
          ! else
          !   do jom = 1, freq%nomeg-1
          !     if ( (freq%freqs(jom) <= w_sc) .and. (w_sc < freq%freqs(jom+1)) ) then
          !       kom = jom
          !       exit
          !     end if
          !   end do
          ! end if
          ! print*, 'iom', iom, w_sc
          ! print*, 'kom', kom, freq%freqs(kom)
          ! Linear interpolation
          ! wdiff = freq%freqs(kom+1)-freq%freqs(kom)
          ! xnm_intp = xnm(kom) + ( xnm(kom+1)-xnm(kom) ) * ( w_sc-freq%freqs(kom) ) / wdiff         

          ! call s1%evaluate( w_sc, 0, f1, iflag)
          ! call s2%evaluate( w_sc, 0, f2, iflag)
          ! xnm_intp = cmplx(f1, f2, 8)

          xnm_intp = reval_aaa_approximant( aaa, cmplx(w_sc, 0.d0, 8) )

          if (ikp==2 .and. iq==2 .and. ie1==1 .and. ie2==1) write(11,'(3f18.6)') w_sc, xnm_intp

          !------------------------
          ! Frequency convolution
          !------------------------

          ! (enk-iu)^2
          zt1 = cmplx( enk, -w_sc, 8)

          sc = zzero
          do jom = 1, freq%nomeg
            zt2 = freq%womeg(jom) / ( freq%freqs(jom)**2 + zt1**2 )
            sc = sc + ( xnm(jom) - xnm_intp ) * zt2
          end do
          sc = sc * zt1 / pi + xnm_intp * sign( 0.5d0, enk )

          ! sum over states
          selfec(ie1,iom,ikp) = selfec(ie1,iom,ikp) + sc

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
end subroutine    
