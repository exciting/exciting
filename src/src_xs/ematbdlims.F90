
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
!!$  select case(typ)
!!$  case(0)
!!$     lo1=1
!!$     hi1=nstsv
!!$     lo2=1
!!$     hi2=nstsv
!!$     n1=nstsv
!!$     n2=nstsv
!!$  case(1)
!!$     lo1=1
!!$     hi1=nstval
!!$     lo2=nstval+1
!!$     hi2=nstsv
!!$     n1=nstval
!!$     n2=nstcon
!!$  case(2)
!!$     lo1=nstval+1
!!$     hi1=nstsv
!!$     lo2=1
!!$     hi2=nstval
!!$     n1=nstcon
!!$     n2=nstval
!!$  case(3)
!!$     lo1=1
!!$     hi1=nstval
!!$     lo2=1
!!$     hi2=nstval
!!$     n1=nstval
!!$     n2=nstval
!!$  case(4)
!!$     lo1=nstval+1
!!$     hi1=nstsv
!!$     lo2=nstval+1
!!$     hi2=nstsv
!!$     n1=nstcon
!!$     n2=nstcon
!!$  end select
!!$  select case(typ)     
!!$  case(0)
!!$     ! all combinations
!!$     lo1=1
!!$     hi1=nstsv
!!$     lo2=1
!!$     hi2=nstsv
!!$     n1=nstsv
!!$     n2=nstsv
!!$  case(1)
!!$     ! 1-2 combinations
!!$     lo1=1
!!$     hi1=istocc0
!!$     lo2=istunocc
!!$     hi2=nstsv
!!$     n1=istocc0-1+1
!!$     n2=nstsv-istunocc+1
!!$  case(2)
!!$     ! 2-1 combinations
!!$     lo1=istunocc0
!!$     hi1=nstsv
!!$     lo2=1
!!$     hi2=istocc
!!$     n1=nstsv-istunocc0+1
!!$     n2=istocc-1+1
!!$  case(3)
!!$     ! 1-1 combinations
!!$     lo1=1
!!$     hi1=istocc0
!!$     lo2=1
!!$     hi2=istocc
!!$     n1=istocc0-1+1
!!$     n2=istocc-1+1
!!$  case(4)
!!$     ! 2-2 combinations
!!$     lo1=istunocc0
!!$     hi1=nstsv
!!$     lo2=istunocc
!!$     hi2=nstsv
!!$     n1=nstsv-istunocc0+1
!!$     n2=nstsv-istunocc+1
!!$  end select
  select case(typ)     
  case(0)
     ! all combinations
     lo1=1
     hi1=nstsv
     lo2=1
     hi2=nstsv
     n1=nstsv
     n2=nstsv
  case(1)
     ! 1-2 combinations
     lo1=1
     hi1=istocc0
     lo2=istunocc0
     hi2=nstsv
     n1=istocc0-1+1
     n2=nstsv-istunocc0+1
  case(2)
     ! 2-1 combinations
     lo1=istunocc0
     hi1=nstsv
     lo2=1
     hi2=istocc0
     n1=nstsv-istunocc0+1
     n2=istocc0-1+1
  case(3)
     ! 1-1 combinations
     lo1=1
     hi1=istocc0
     lo2=1
     hi2=istocc0
     n1=istocc0-1+1
     n2=istocc0-1+1
  case(4)
     ! 2-2 combinations
     lo1=istunocc0
     hi1=nstsv
     lo2=istunocc0
     hi2=nstsv
     n1=nstsv-istunocc0+1
     n2=nstsv-istunocc0+1
  end select
end subroutine ematbdlims
