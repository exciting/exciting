! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created by Ronaldo Rodrigues Pela, May 2019
! Reference: https://arxiv.org/abs/2102.02630

!> Module that manages what concerns charge density in RT-TDDFT calculations
module rttddft_Density
  implicit none

  private

  public :: UpdateDensity

contains
  !> In UpdateDensity, we obtain the charge density at time t
  !> It is calculated from the WFs, using the same scheme as in the GS
  !> calcultations (refer to scf_cycle.f90 for the case
  !> input%groundstate%useDensityMatrix .false.)
  !> @param[in]   it            number of the current iteration
  !>                            (employed to give possible warnings)
  !> @param[in]   timeini       time (in seconds) elapsed since exciting was started
  !>                            (employed for timing)
  !> @param[out]   timefinal    time (in seconds) after executing this subroutine
  !> @param[out]   timerho      time (in seconds) taken to execute rhovalk,
  !>                            genrhoir, and eventually mpisumrhoandmag
  !> @param[out]   timesymrf    time (in seconds) taken to execute symrf
  !> @param[out]   timerfmtctof time (in seconds) taken to execute rfmtctof
  !> @param[out]   timerhocr    time (in seconds) taken to obtain and add the
  !>                            density of core electrons
  !> @param[out]   timecharge   time (in seconds) taken to execute charge
  !> @param[out]   timerhonorm  time (in seconds) taken to execute rhonorm
  subroutine UpdateDensity( it, timeini, timefinal, timerho, &
      & timesymrf, timerfmtctof, timerhocr, timecharge, timerhonorm )
    use modmpi
    use precision, only: dp
    use modmain, only : iscl
    use modinput, only: input
    use mod_kpoint, only: nkpt
    use mod_potential_and_density, only: rhomt, rhoir
    use rttddft_GlobalVariables, only: evecfv_time, evecsv

    implicit none

    !> it:  number of the current iteration (employed to give possible warnings)
    integer, intent(in)             :: it
    !> timeini: time (in seconds) elapsed since exciting was started (employed for timing)
    real(dp),intent(in), optional   :: timeini
    !> timefinal: time (in seconds) after executing this subroutine
    real(dp),intent(out), optional  :: timefinal
    !> timerho: time (in seconds) taken to execute rhovalk, genrhoir, and eventually mpisumrhoandmag
    !> timesymrf: time (in seconds) taken to execute symrf
    !> timerfmtctof: time (in seconds) taken to execute rfmtctof
    !> timerhocr:  time (in seconds) taken to obtain and add the density of core electrons
    real(dp),intent(out), optional  :: timerho, timesymrf, timerfmtctof, timerhocr
    !> timecharge:   time (in seconds) taken to execute charge
    !> timerhonorm:  time (in seconds) taken to execute rhonorm
    real(dp),intent(out), optional  :: timecharge, timerhonorm

    integer                         :: ik, first_kpt, last_kpt
    real(dp)                        :: timei,timef


    if( present( timeini ) ) timei = timeini
    rhomt(:,:,:) = 0._dp
    rhoir(:) = 0._dp

#ifdef MPI
  first_kpt = firstk(rank)
  last_kpt = lastk(rank)
#else
  first_kpt = 1
  last_kpt = nkpt
#endif

    do ik = first_kpt, last_kpt
      call rhovalk( ik, evecfv_time(:, :, ik), evecsv(:, :, ik) )
      call genrhoir( ik, evecfv_time(:, :, ik), evecsv(:, :, ik) )
    end do
#ifdef MPI
    call mpisumrhoandmag
#endif
    if ( present( timerho ) ) then
      call timesec( timef )
      timerho = timef - timei
      timei = timef
    end if
    ! symmetrise the density
    call symrf( input%groundstate%lradstep, rhomt, rhoir )
    if( present( timesymrf ) ) then
      call timesec( timef )
      timesymrf = timef - timei
      timei = timef
    end if
    ! convert the density from a coarse to a fine radial mesh
    call rfmtctof(rhomt)
    if(present(timerfmtctof)) then
      call timesec(timef)
      timerfmtctof = timef-timei
      timei = timef
    end if
    ! generate the core wavefunctions and densities
    !call gencore
    ! add the core density to the total density
    call addrhocr
    if( present( timerhocr ) ) then
      call timesec( timef )
      timerhocr = timef - timei
      timei = timef
    end if
    ! calculate the charges
    iscl = it
    call charge
    if(present(timecharge)) then
      call timesec(timef)
      timecharge = timef-timei
      timei = timef
    end if
    ! normalise the density
    if ( input%xs%rt_tddft%normalize_WF ) then
      call rhonorm
    end if
    if(present(timefinal)) then
      call timesec(timef)
      timefinal = timef
    end if
    if(input%xs%rt_tddft%normalize_WF .and. present(timerhonorm)) then
      timerhonorm = timef-timei
    end if
  end subroutine updatedensity

end module rttddft_Density
