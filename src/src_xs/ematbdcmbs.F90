! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine ematbdcmbs(etyp)
  use modxs, only: nst1, nst2, nst3, nst4,&
                & istl1, istl2, istl3, istl4,&
                & istu1, istu2, istu3, istu4

  implicit none

  ! Arguments
  integer, intent(in) :: etyp
  
  select case(etyp)
    case(0)
      ! All combinations
      ! nst1=nstsv, istl1=1, istu1=nstsv
      ! nst2=nstsv, istl2=1, istu2=nstsv
      call ematbdlims(0, nst1, istl1, istu1, nst2, istl2, istu2)
      nst3 = 0
      nst4 = 0
    case(1)
      ! o-u combinations
      ! nst1=sto1-sta1+1, istl1=sta1, istu1=sto1,
      ! nst2=sto2-sta2+1, istl2=istunocc0+sta2-1, istu2=istunocc0+sto2-1
      call ematbdlims(1, nst1, istl1, istu1, nst2, istl2, istu2)
      ! u-o combinations
      ! nst3=sto2-sta2+1, istl3=istunocc0+sta2-1, istu3=istunocc0+sto2-1
      ! nst4=sto1-sta1+1, istl4=sta1, istu4=sto1,
      call ematbdlims(2, nst3, istl3, istu3, nst4, istl4, istu4)
    case(2)
      ! o-o combinations
      ! nst1=sto1-sta1+1, istl1=sta1, istu1=sto1,
      ! nst2=sto1-sta1+1, istl2=sta1, istu2=sto1
      call ematbdlims(3, nst1, istl1, istu1, nst2, istl2, istu2)
      ! u-u combinations
      ! nst3=sto2-sta2+1, istl3=istunocc0+sta2-1, istu3=istunocc0+sto2-1
      ! nst4=sto2-sta2+1, istl4=istunocc0+sta2-1, istu4=istunocc0+sto2-1
      call ematbdlims(4, nst3, istl3, istu3, nst4, istl4, istu4)
  end select

end subroutine ematbdcmbs
