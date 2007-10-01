
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrvfir(sym,rvfir,srvfir)
use modmain
implicit none
! arguments
integer, intent(in) :: sym(3,3)
real(8), intent(in) :: rvfir(ngrtot,3)
real(8), intent(out) :: srvfir(ngrtot,3)
! local variables
integer i,ir
real(8) s(3,3),v(3)
! convert symmetry matrix from lattice to Cartesian coordinates
s(:,:)=dble(sym(:,:))
call r3mm(s,ainv,s)
call r3mm(avec,s,s)
! perform parallel transport rotation on the vector function
do i=1,3
  call symrfir(sym,rvfir(1,i),srvfir(1,i))
end do
! rotate the vectors at each point
do ir=1,ngrtot
  v(:)=srvfir(ir,:)
  call r3mv(s,v,v)
  srvfir(ir,:)=v(:)
end do
return
end subroutine
