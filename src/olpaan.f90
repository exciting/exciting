
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaan(is,ia,ngp,apwalm)
use modmain
use modfvsystem
implicit none
! arguments

integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)


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
    if(packed) then
    
      call ZHPR2 ( 'U', ngp, zhalf, conjg(apwalm(1,io,lm,ias)), 1, conjg(apwalm(1,io,lm,ias)), 1, op )
    else
      call ZHER2 ( 'U', ngp, zhalf, conjg(apwalm(1,io,lm,ias)), 1, conjg(apwalm(1,io,lm,ias)), 1, o,ohrank)
    endif
    end do
  end do
end do
return
end subroutine

