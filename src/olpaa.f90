
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaa(tapp,is,ia,ngp,apwalm,v,o)
use modmain
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: v(nmatmax)
complex(8), intent(inout) :: o(*)
! local variables
integer ias,l,m,lm,io
! external functions
complex(8) zdotu
external zdotu
ias=idxas(ia,is)
do l=0,lmaxmat
  do m=-l,l
    lm=idxlm(l,m)
    do io=1,apword(l,is)
      call zmatinp(tapp,ngp,zhalf,apwalm(1,io,lm,ias),apwalm(1,io,lm,ias),v,o)
    end do
  end do
end do
return
end subroutine

