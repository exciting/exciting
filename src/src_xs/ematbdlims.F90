! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine ematbdlims(typ, n1, lo1, hi1, n2, lo2, hi2)
  use mod_misc, only: task
  use mod_eigenvalue_occupancy, only: nstsv
  use modxs, only: sta1, sto1, sta2, sto2, istunocc0, istocc0

  implicit none

  ! Arguments
  integer, intent(in) :: typ
  integer, intent(out) :: n1, n2, lo1, hi1, lo2, hi2

  if(task .ge. 440 .and. task .le. 446) then

    select case(typ)
      case(0)
        ! All combinations
        lo1 = 1
        hi1 = nstsv
        lo2 = 1
        hi2 = nstsv
        n1 = nstsv
        n2 = nstsv
      case(1)
        ! o-u combinations
        lo1 = sta1
        hi1 = sto1
        lo2 = istunocc0+sta2-1
        hi2 = istunocc0+sto2-1
        n1 = sto1- sta1 + 1
        n2 = sto2 -sta2 + 1
      case(2)
        ! u-o combinations
        lo1 = istunocc0+sta2-1
        hi1 = istunocc0+sto2-1
        lo2 = sta1
        hi2 = sto1
        n1 = sto2-sta2+1
        n2 = sto1-sta1+1
      case(3)
        ! o-o combinations
        lo1 = sta1
        hi1 = sto1
        lo2 = sta1
        hi2 = sto1
        n1 = sto1-sta1+1
        n2 = sto1-sta1+1
      case(4)
        ! u-u combinations
        lo1 = istunocc0+sta2-1
        hi1 = istunocc0+sto2-1
        lo2 = istunocc0+sta2-1
        hi2 = istunocc0+sto2-1
        n1 = sto2-sta2+1
        n2 = sto2-sta2+1
    end select

  else

    select case(typ)
      case(0)
        ! All combinations
        lo1 = 1
        hi1 = nstsv
        lo2 = 1
        hi2 = nstsv
        n1 = nstsv
        n2 = nstsv
      case(1)
        ! o-u combinations
        lo1 = 1
        hi1 = istocc0
        lo2 = istunocc0
        hi2 = nstsv
        n1 = istocc0 - 1 + 1
        n2 = nstsv - istunocc0 + 1
      case(2)
        ! u-o combinations
        lo1 = istunocc0
        hi1 = nstsv
        lo2 = 1
        hi2 = istocc0
        n1 = nstsv - istunocc0 + 1
        n2 = istocc0 - 1 + 1
      case(3)
        ! o-o combinations
        lo1 = 1
        hi1 = istocc0
        lo2 = 1
        hi2 = istocc0
        n1 = istocc0 - 1 + 1
        n2 = istocc0 - 1 + 1
      case(4)
        ! u-u combinations
        lo1 = istunocc0
        hi1 = nstsv
        lo2 = istunocc0
        hi2 = nstsv
        n1 = nstsv - istunocc0 + 1
        n2 = nstsv - istunocc0 + 1
    end select

  end if

end subroutine ematbdlims
