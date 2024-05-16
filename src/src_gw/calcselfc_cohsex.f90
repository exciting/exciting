!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! in static COHSEX model
!==================================================================

subroutine calcselfc_cohsex(ikp,iq,mdim)
    use modinput, only: input
    use mod_eigenvalue_occupancy, only: nstfv
    use modgw,   only : ibgw, nbgw, fdebug
    use mod_bands, only: nomax
    use mod_selfenergy, only: selfec, sigsx, sigch, mwm
#include "mod_gw_degeneracies.inc"
    use mod_gw_degeneracies, only: get_degenerate_limits_qp_interval_ikp, &
                                   degenerate_subspaces
    use precision, only: i32, dp
    use constants, only: zzero

    ! input variables
    implicit none
    integer(i32), intent(in) :: ikp
    integer(i32), intent(in) :: iq
    integer(i32), intent(in) :: mdim
    ! local variables            
    integer(i32) :: ie1, ie2
    complex(dp)  :: sc

    ! For the averaging over degenerate states
    integer(i32) :: ispace_init, ispace_final, ispace, lowband, upband, size_deg

    ! First obtain the limits for degenerate subspaces for the given irreducible point. 
    call get_degenerate_limits_qp_interval_ikp(ikp, ispace_init, ispace_final)

    ! Now iterate over the degenerate spaces
    do ispace = ispace_init, ispace_final 

      lowband  = degenerate_subspaces(1, ispace, ikp)
      upband   = degenerate_subspaces(2, ispace, ikp)
      size_deg = degenerate_subspaces(3, ispace, ikp)

      do ie1 = lowband, upband
        sc = zzero
        do ie2 = 1, mdim
          ! occupied states
          if ((ie2<=nomax).or.(ie2>nstfv)) then
            sc = sc - 0.5_dp*mwm(ie1,ie2,1)
            ! \Sigma_{SEX}
            sigsx(QP_ADJUST_RANGE(lowband,upband),ikp) = & 
              sigsx(QP_ADJUST_RANGE(lowband,upband),ikp) - mwm(ie1,ie2,1) / size_deg
          else
            sc = sc + 0.5_dp*mwm(ie1,ie2,1)
          end if
          ! \Sigma_{COH}
          sigch(QP_ADJUST_RANGE(lowband,upband),ikp) = & 
            sigch(QP_ADJUST_RANGE(lowband,upband),ikp) + 0.5_dp*mwm(ie1,ie2,1) / size_deg
        end do ! ie2
        selfec(QP_ADJUST_RANGE(lowband,upband),1,ikp) = & 
          selfec(QP_ADJUST_RANGE(lowband,upband),ikp,1) + sc / size_deg
      end do ! ie1
    end do
    
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
