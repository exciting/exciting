! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created by Ronaldo Rodrigues Pela, May 2019
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module that manages what concerns charge density in RT-TDDFT calculations
module rttddft_Density
  use asserts, only: assert
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_density, timesec_RTTDDFT

  implicit none

  private

  public :: UpdateDensity

contains
  !> In UpdateDensity, we obtain the charge density at time t
  !> It is calculated from the WFs, using the same scheme as in the GS
  !> calcultations (refer to scf_cycle.f90 for the case
  !> input%groundstate%useDensityMatrix .false.)
  subroutine UpdateDensity( first_kpt, wavefunctions, wavefunctions_sv, &
    it, normalize, l_rad_step, printTimings, t_dens )
    use modmpi, only: mpi_env_k, distribute_loop
    use precision, only: dp, i32
    use modmain, only : iscl
    use mod_kpoint, only: nkpt
    use mod_potential_and_density, only: rhomt, rhoir

    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> Wavefunctions in LAPW basis (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), intent(in) :: wavefunctions(:, :, first_kpt :)
    !> '2nd variation' wavefunctions (nstsv, nstsv, first_kpt : last_kpt)
    complex(dp), intent(in) :: wavefunctions_sv(:, :, first_kpt :)
    !> number of the current iteration (employed to give possible warnings)
    integer, intent(in)             :: it
    !> If `.true.`, normalize the charge density
    logical, intent(in)             :: normalize
    !> radial step length
    integer(i32), intent(in)        :: l_rad_step
    !> Object that packs information about printing of timings [[Print_Timings]]
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the electronic density
    type(Timing_RTTDDFT_density), optional, intent(out) :: t_dens

    integer(i32)                    :: ik, last_kpt
    real(dp)                        :: ti, tstart 
    logical                         :: timings_general, timings_detailed

    timings_general = .false.
    timings_detailed = .false.
    if( present( printTimings ) ) then
      call assert( present( t_dens ), 'Optional argument t_dens must also be present' )
      call printTimings%get( timings_general, timings_detailed )
      if( timings_detailed ) call assert( timings_general, &
        'timings_general must be true if timings_detailed is true' )
    end if
    if( timings_general ) then
      call timesec( ti )
      tstart = ti
    end if

    rhomt(:, :, :) = 0._dp
    rhoir(:) = 0._dp
    last_kpt = ubound( wavefunctions, 3 )
    
    ! rhovalk has omp critical inside, and we use reduction for rhoir
#ifdef USEOMP
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik) &
    !$OMP SHARED(first_kpt, last_kpt, wavefunctions, wavefunctions_sv)
    !$OMP DO
#endif
    do ik = first_kpt, last_kpt
      call rhovalk( ik, wavefunctions(:, :, ik), wavefunctions_sv(:, :, ik) )
      call genrhoir( ik, wavefunctions(:, :, ik), wavefunctions_sv(:, :, ik) )
    end do
#ifdef USEOMP
    !$OMP END PARALLEL
#endif
    
#ifdef MPI
    call mpisumrhoandmag( mpi_env_k )
#endif

    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%rho )

    ! symmetrise the density
    call symrf( l_rad_step, rhomt, rhoir )
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%symrf )

    ! convert the density from a coarse to a fine radial mesh
    call rfmtctof( rhomt )
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%rfmtctof )

    ! generate the core wavefunctions and densities
    !call gencore
    ! add the core density to the total density
    call addrhocr()
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%addrhocr )

    ! calculate the charges
    iscl = it
    call charge()
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%charge )

    ! normalise the density
    if ( normalize ) then
      call rhonorm()
      if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%rhonorm )
    end if
    
    if( timings_general ) call timesec_RTTDDFT( tstart, t_dens%total )

  end subroutine updatedensity

end module rttddft_Density
