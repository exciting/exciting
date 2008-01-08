
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematbdlims(typ,n1,lo1,hi1,n2,lo2,hi2)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: typ
  integer, intent(out) :: n1,n2,lo1,hi1,lo2,hi2
  ! set limits for band combinatios
  select case(typ)
  case(0)
     lo1=1
     hi1=nstsv
     lo2=1
     hi2=nstsv
     n1=nstsv
     n2=nstsv
  case(1)
     lo1=1
     hi1=nstval
     lo2=nstval+1
     hi2=nstsv
     n1=nstval
     n2=nstcon
  case(2)
     lo1=nstval+1
     hi1=nstsv
     lo2=1
     hi2=nstval
     n1=nstcon
     n2=nstval
  case(3)
     lo1=1
     hi1=nstval
     lo2=1
     hi2=nstval
     n1=nstval
     n2=nstval
  case(4)
     lo1=nstval+1
     hi1=nstsv
     lo2=nstval+1
     hi2=nstsv
     n1=nstcon
     n2=nstcon
  end select
end subroutine ematbdlims
