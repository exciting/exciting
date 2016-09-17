!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! in static COHSEX model
!==================================================================
subroutine calcselfc_cohsex(ikp,iq,mdim)
    use modinput
    use modmain, only : zzero, nstsv
    use modgw,   only : ibgw, nbgw, nomax, selfec, mwm, sigch, sigsx, &
    &                   fdebug
    ! input variables
    implicit none
    integer(4), intent(in) :: ikp
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: mdim
    ! local variables            
    integer(4) :: ie1, ie2
    complex(8) :: sc
    
    do ie1 = ibgw, nbgw
      sc = zzero
      do ie2 = 1, mdim
        ! occupied states
        if ((ie2<=nomax).or.(ie2>nstsv)) then
          sc = sc-0.5d0*mwm(ie1,ie2,1)
          ! \Sigma_{SEX}
          sigsx(ie1,ikp) = sigsx(ie1,ikp)-mwm(ie1,ie2,1)
        else
          sc = sc+0.5d0*mwm(ie1,ie2,1)
        end if
        ! \Sigma_{COH}
        sigch(ie1,ikp) = sigch(ie1,ikp)+0.5d0*mwm(ie1,ie2,1)
      end do ! ie2
      selfec(ie1,1,ikp) = selfec(ie1,ikp,1)+sc
    end do ! ie1
    
    if (input%gw%debug) then
      write(fdebug,*) 'COHSEX: CORRELATION SELF-ENERGY: iq=', iq, ' ikp=', ikp
      write(fdebug,*) 'state    Sigma_c'
      do ie1 = ibgw, nbgw
          write(fdebug,*) ie1, selfec(ie1,1,ikp)
      end do
      write(fdebug,*)
    end if

    return 
end subroutine
