! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Created April 2003 (JKD)

!> Computes the effective potential by adding together the Coulomb and
!> exchange-correlation potentials. See routines [[potcoul]] and [[potxc]].
subroutine poteff( calc_xc )
  use constants, only: y00
  use mod_atoms, only: nspecies, natoms, idxas
  use mod_muffin_tin, only: nrmt, lmmaxvr, lmmaxinr, nrmtinr
  use mod_potential_and_density, only: veffmt, veffir, vhalfmt, vhalfir, vclmt, &
    vclir, vxcmt, vxcir
  use mod_timing, only: stopwatch, timepot
  use modinput, only: input
  use precision, only: i32, dp

  implicit none
  !> If `.true.`, XC part of the effective potential should be obtained
  logical, intent(in) :: calc_xc

  integer(i32) :: is, ia, ias, ir, lmmax
  real(dp) :: shift, ts0, ts1
  logical :: dfthalf_on

  call stopwatch( "exciting:poteff", 1 ) 
  call timesec ( ts0 )

  shift = input%groundstate%energyref
  dfthalf_on = associated( input%groundstate%dfthalf )
  if( dfthalf_on ) dfthalf_on = .not. input%groundstate%dfthalf%NSCF

  veffmt = 0._dp
  veffir = 0._dp
  
  ! compute the exchange-correlation potential
  if ( calc_xc ) call potxc()
  ! compute the Coulomb potential
  call potcoul()

  ! add Coulomb and exchange-correlation potentials together
  ! muffin-tin part
  vclmt(1, :, :) = vclmt(1, :, :) + shift / y00
  do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas (ia, is)
      lmmax = lmmaxinr
      do ir = 1, nrmt (is)
        if ( ir > nrmtinr(is) ) lmmax = lmmaxvr
        veffmt(1 : lmmax, ir, ias) = vclmt(1 : lmmax, ir, ias) + vxcmt(1 : lmmax, ir, ias)
        if ( dfthalf_on ) veffmt(1 : lmmax, ir, ias) = veffmt(1 : lmmax, ir, ias) + vhalfmt (1 : lmmax, ir, ias)
      end do
    end do
  end do
  ! interstitial part
  veffir = vclir + vxcir + shift
  if ( dfthalf_on ) veffir = veffir + vhalfir
  
  call timesec ( ts1 )
  timepot = timepot + ts1 - ts0
  call stopwatch( "exciting:poteff", 0)
end subroutine
