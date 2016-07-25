!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! using the frequency convolution
!==================================================================
subroutine calcselfc_freqconv(ikp,iq,mdim)
    use modinput
    use modmain, only : pi, zzero, evalsv, idxas, evalcr, efermi
    use modgw,   only : ibgw, nbgw, nstse, kset, kqset, freq, selfec, mwm, &
    &                   ncg, corind, fdebug
    ! input variables
    implicit none
    integer(4), intent(in) :: ikp
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: mdim
    ! local variables            
    integer(4) :: ik, jk, jkp
    integer(4) :: ia, is, ias, ic, icg
    integer(4) :: ie1, ie2
    integer(4) :: iom, jom
    real(8)    :: enk
    complex(8) :: xnm(1:freq%nomeg)
    complex(8) :: sc, zt1, zt2, zt3
    
    ! k point
    ik = kset%ikp2ik(ikp)
    ! k-q point
    jk = kqset%kqid(ik,iq)
    jkp = kset%ik2ikp(jk)

    !------------------------
    ! loop over frequencies
    !------------------------
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iom,ie1,ie2,enk,xnm,icg,is,ia,ic,ias,zt1,sc,jom,zt2,zt3)
!$OMP DO
#endif
    do iom = 1, freq%nomeg
      ! loop over states
      do ie1 = ibgw, nbgw
      
        ! sum over states
        do ie2 = 1, mdim
        
          if (ie2<=nstse) then
            !============================= 
            ! Valence electron contribution
            !============================= 
            enk = evalsv(ie2,jkp)-efermi
            xnm(1:freq%nomeg) = mwm(ie1,ie2,1:freq%nomeg)
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
            xnm(1:freq%nomeg) = mwm(ie1,ie2,1:freq%nomeg)
          end if ! val/cor
          
          !------------------------
          ! Frequency convolution
          !------------------------
          ! (enk-iu)^2
          zt1 = cmplx(enk,-freq%freqs(iom),8)
          sc = 0.d0
          do jom = 1, freq%nomeg
            zt2 = cmplx(freq%freqs(jom)*freq%freqs(jom),0.0d0,8)
            zt3 = freq%womeg(jom)/(zt1*zt1+zt2)
            sc = sc+(xnm(jom)-xnm(iom))*zt3
          end do
          sc = zt1*sc/pi+xnm(iom)*sign(0.5d0,enk)

          ! sum over states
          selfec(ie1,iom,ikp) = selfec(ie1,iom,ikp)+sc
                  
        end do ! ie2
        
      end do ! ie1
    end do ! iom
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    

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
