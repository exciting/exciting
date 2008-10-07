
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

logical function transijst(ik,ist,jst,trans)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: ik,ist,jst,trans(3,ndftrans)
  ! local variables
  integer :: l,ikt,it,jt,ir,jr
  transijst=.false.
  do l=1,ndftrans
     ikt=dftrans(1,l)
     it=dftrans(2,l)
     jt=dftrans(3,l)
     if ((ikt.eq.0).or.(ikt.eq.ik)) then
        ! it-jt, it-x, x-jt and x-x combinations
        if (((it.eq.0).or.(it.eq.ist)).and.((jt.eq.0).or.(jt.eq.jst))) then
           transijst=.true.
           return
        end if
        ! ranges between |it| and |jt| (for negative it and jt)
        ir=-it
        jr=-jt
        if (((ir.gt.0).and.(jr.gt.0)).and.(ir.lt.jr)) then
           if ((ist.ge.ir).and.(jst.le.jr)) transijst=.true.
           return
        end if
     end if
  end do
end function transijst
