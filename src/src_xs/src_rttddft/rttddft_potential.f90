!> Module that manages what concerns the KS potential in RT-TDDFT calculations
module rttddft_potential
#include "asserts.fpp"
  use precision, only: dp
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_potential, timesec_RTTDDFT

  implicit none

  private

  public :: update_potential

contains
  !> Obtain the KS potential given the electron density 
  subroutine update_potential( printTimings, t_pot, coulomb_only )
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the KS potential
    type(Timing_RTTDDFT_potential), optional, intent(inout) :: t_pot
    !> Whether effective potential without the XC part should be calculated
    logical, optional, intent(in) :: coulomb_only
    
    logical :: timings_general, timings_detailed, calc_xc
    real(dp) :: ti, t_start

    timings_general = .false.
    timings_detailed = .false.
    if( present( printTimings ) ) then 
      CALL_ASSERT( present( t_pot ), 't_pot must be present if printTimings is')
      call printTimings%get( timings_general, timings_detailed )
    end if 
    if( timings_general ) then
      call timesec( ti )
      t_start = ti
    end if

    calc_xc = .true.
    if ( present( coulomb_only ) ) calc_xc = ( .not. coulomb_only )

    ! Compute the effective potential (with the updated density)
    call poteff( calc_xc )

    if( timings_detailed ) call timesec_RTTDDFT( ti, t_pot%poteff )
    ! Fourier transform effective potential to G-space
    call genveffig()
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_pot%genveffig )
    call genmeffig()
    if( timings_general ) then
      call timesec_RTTDDFT( t_start, t_pot%total )
      if( timings_detailed ) t_pot%genmeffig = t_start - ti
    end if
  end subroutine

end module