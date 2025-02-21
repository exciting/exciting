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
  use precision, only: dp, i32
  use modmpi, only: mpi_env_k
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_density, timesec_RTTDDFT
  use constants, only: zzero, real_zero
  use mod_convergence, only : iscl
  use mod_potential_and_density, only: rhomt, rhoir
  use mod_eigenvalue_occupancy, only: occsv, nstfv
  use rttddft_Wavefunction, only: wavefunction_set
  use mod_rhoir, only: genrhoir
  use mod_rhovalk, only: rhovalk

  implicit none

  private

  public :: update_density, save_and_frozen, frozen, groundstate

  !> Enum with the density evaluation mode
  !> There are 4 options: get active density, save density, ground density
  enum, bind(C)
    enumerator :: density_case
    enumerator :: active_and_frozen, save_and_frozen, groundstate, frozen
  end enum

contains
  !> In `update_density`, we obtain the charge density at time \(t\) 
  !> and put it in global 'rho_mt" and 'rho_ir' arrays.
  !> It is calculated from the WFs, using the same scheme as in the GS
  !> calcultations (refer to `scf_cycle.f90` for the case
  !> `input%groundstate%useDensityMatrix` `.false.`)
  subroutine update_density( first_kpt, psi, it, normalize, l_rad_step, &
      rhomt_frozen, rhoir_frozen, printTimings, t_dens, dens_case )
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> Set of KS wavefunctions
    class(wavefunction_set), intent(in) :: psi
    !> number of the current iteration (employed to give possible warnings)
    integer(i32), intent(in) :: it
    !> If `.true.`, normalize the charge density
    logical, intent(in) :: normalize
    !> radial step length
    integer(i32), intent(in) :: l_rad_step
    !> Frozen part of the muffin-tin density (lmmaxvr, nrmtmax, natmtot)
    real(dp), optional, intent(in) :: rhomt_frozen(:, :, :)
    !> Frozen part of the IR density (ngrtot)
    real(dp), optional, intent(in) :: rhoir_frozen(:)
    !> Object that packs information about printing of timings [[Print_Timings]]
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the electronic density
    type(Timing_RTTDDFT_density), optional, intent(out) :: t_dens
    !> Enum telling which part of the wavefuntion should be used for density evaluation
    integer(kind( density_case )), optional, intent(in) :: dens_case

    integer(i32) :: ik, i, last_kpt, first_active
    real(dp) :: ti, tstart 
    logical :: timings_general, timings_detailed, add_frozen
    integer(kind( density_case )) :: dens_case_
    real(dp), allocatable :: occupations(:, :)

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

    rhomt = real_zero
    rhoir = real_zero

    first_active = psi%first_active()

    dens_case_ = active_and_frozen
    if ( present( dens_case ) ) dens_case_ = dens_case
    add_frozen = .false.
    
    if ( present( rhomt_frozen ) .or. present( rhoir_frozen ) ) then
      call assert( present( rhomt_frozen ) .and. present( rhoir_frozen ), &
        'Both contributions to frozen density should be provided to update_density' )
      ! for the ground state it is convenient to ignore precalculated frozen density
      add_frozen = ( .not. ( dens_case_ == groundstate ) )
    end if

    if ( dens_case_ == frozen ) then
      call assert( .not. add_frozen, &
        'frozen density requested with frozen rhoir and rhomt present in update_density' )
      call assert( psi%has_frozen(), &
        'frozen density requested with no frozen wavefunctions' )
    end if

    last_kpt = first_kpt + psi%n_kpts() - 1
    select case ( dens_case_ )
    case( active_and_frozen, save_and_frozen )
      allocate( occupations, source = occsv(first_active : nstfv, first_kpt : last_kpt))
    case( frozen )
      allocate( occupations, source = occsv(1 : first_active - 1, first_kpt : last_kpt))
    case( groundstate )
      allocate( occupations, source = occsv(:, first_kpt : last_kpt))
    case default
      call assert( .false., 'unknown dens_case_' )
    end select

    ! rhovalk and rhoir have omp critical inside
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, ik) &
    !$OMP SHARED(first_kpt, psi, occupations, dens_case_, rhomt, rhoir)
    !$OMP DO
    do i = 1, size( occupations, 2 )
      ik = first_kpt + i - 1
      select case ( dens_case_ )
      case( active_and_frozen )
        call rhovalk( ik, psi%active(:, :, i), occupations(:, i), rhomt )
        call genrhoir( ik, psi%active(:, :, i), occupations(:, i), rhoir )
      case( save_and_frozen )
        call rhovalk( ik, psi%active_save(:, :, i), occupations(:, i), rhomt )
        call genrhoir( ik, psi%active_save(:, :, i), occupations(:, i), rhoir )
      case( frozen )
        call rhovalk( ik, psi%frozen(:, :, i), occupations(:, i), rhomt )
        call genrhoir( ik, psi%frozen(:, :, i), occupations(:, i), rhoir )
      case( groundstate )
        call rhovalk( ik, psi%groundstate(:, :, i), occupations(:, i), rhomt )
        call genrhoir( ik, psi%groundstate(:, :, i), occupations(:, i), rhoir )
      end select
    end do
    !$OMP END PARALLEL

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

    if ( add_frozen ) then
      ! frozen density arrays contain core contribution along with the valence one
      rhoir = rhoir + rhoir_frozen
      rhomt = rhomt + rhomt_frozen
    else
      ! add the core density to the total density
      call addrhocr()
      if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%addrhocr )
    end if

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

  end subroutine update_density

end module rttddft_Density
