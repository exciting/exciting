

! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine olpaan(overlap, is, ia, ngp, apwalm)
use modmain
use modinput
use modfvsystem
implicit none
! arguments
type (hermiteanmatrix), intent(inout)::overlap
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax, apwordmax, lmmaxapw, natmtot)
complex(8)::x(ngp), y(ngp)

! local variables
integer::ias, l, m, lm, io
! external functions
complex(8) zdotu
external zdotu
ias=idxas(ia, is)
do l=0, input%groundstate%lmaxmat
  do m=-l, l
    lm=idxlm(l, m)
    do io=1, apword(l, is)
    x=conjg(apwalm(1:ngp, io, lm, ias))
    y=conjg(apwalm(1:ngp, io, lm, ias))
    call  Hermiteanmatrix_rank2update(overlap, ngp, zhalf, x, y)
    end do
  end do
end do
return
end subroutine
