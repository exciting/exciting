! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created by Ronaldo Rodrigues Pela, May 2019
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

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
  subroutine UpdateDensity( it, timeini, timefinal, timerho, &
      & timesymrf, timerfmtctof, timerhocr, timecharge, timerhonorm )
    use modmpi, only: mpi_env_k, distribute_loop
    use precision, only: dp
    use modmain, only : iscl
    use modinput, only: input
    use mod_kpoint, only: nkpt
    use mod_potential_and_density, only: rhomt, rhoir
    use rttddft_GlobalVariables, only: evecfv_time, evecsv

    implicit none

    !> number of the current iteration (employed to give possible warnings)
    integer, intent(in)             :: it
    !> time (in seconds) elapsed since exciting was started (employed for timing)
    real(dp),intent(in), optional   :: timeini
    !> time (in seconds) after executing this subroutine
    real(dp),intent(out), optional  :: timefinal
    !> time (in seconds) taken to execute `rhovalk` (which generates the 
    !> valence charge density from the eigenvectors for a certain `k-point` 
    !> inside the Muffin-Tins),`genrhoir` (which does the same as `rhovalk`, but
    !> for the interstitial region), and eventually `mpisumrhoandmag` (which 
    !> sums the charge density of all MPI processes)
    real(dp),intent(out), optional  :: timerho
    !> time (in seconds) taken to execute `symrf`
    real(dp),intent(out), optional  :: timesymrf
    !> time (in seconds) taken to execute `rfmtctof`
    real(dp),intent(out), optional  :: timerfmtctof
    !> time (in seconds) taken to obtain and add the density of core electrons
    real(dp),intent(out), optional  :: timerhocr
    !> time (in seconds) taken to execute `charge` 
    real(dp),intent(out), optional  :: timecharge
    !> time (in seconds) taken to execute `rhonorm`
    real(dp),intent(out), optional  :: timerhonorm

    integer                         :: ik, first_kpt, last_kpt
    real(dp)                        :: timei,timef


    if( present( timeini ) ) timei = timeini
    rhomt(:,:,:) = 0._dp
    rhoir(:) = 0._dp

    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

    do ik = first_kpt, last_kpt
      call rhovalk( ik, evecfv_time(:, :, ik), evecsv(:, :, ik) )
      call genrhoir( ik, evecfv_time(:, :, ik), evecsv(:, :, ik) )
    end do
    
#ifdef MPI
    call mpisumrhoandmag(mpi_env_k)
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
    if ( input%xs%realTimeTDDFT%normalizeWF ) then
      call rhonorm
    end if
    if(present(timefinal)) then
      call timesec(timef)
      timefinal = timef
    end if
    if(input%xs%realTimeTDDFT%normalizeWF .and. present(timerhonorm)) then
      timerhonorm = timef-timei
    end if
  end subroutine updatedensity

end module rttddft_Density
