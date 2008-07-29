
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematbdcmbs(etyp)
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: etyp
  select case(etyp)
  case(0)
     ! all combinations
     call ematbdlims(0,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
     nst3=0; nst4=0
  case(1)
     ! o-u combinations
     call ematbdlims(1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
     ! u-o combinations
     call ematbdlims(2,nst3,istlo3,isthi3,nst4,istlo4,isthi4)
  case(2)
     ! o-o combinations
     call ematbdlims(3,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
     ! u-u combinations
     call ematbdlims(4,nst3,istlo3,isthi3,nst4,istlo4,isthi4)
  end select
end subroutine ematbdcmbs
