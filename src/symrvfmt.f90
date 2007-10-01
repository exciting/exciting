
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrvfmt(lrstp,is,sym,rvfmt,srvfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: lrstp
integer, intent(in) :: is
integer, intent(in) :: sym(3,3)
real(8), intent(in) :: rvfmt(lmmaxvr,nrmtmax,3)
real(8), intent(out) :: srvfmt(lmmaxvr,nrmtmax,3)
! local variables
integer i,ir,lmmax,lm
real(8) s(3,3),v(3)
! perform parallel transport rotation on the vector function
do i=1,3
  call symrfmt(lrstp,is,sym,rvfmt(1,1,i),srvfmt(1,1,i))
end do
! convert symmetry matrix from lattice to Cartesian coordinates
s(:,:)=dble(sym(:,:))
call r3mm(s,ainv,s)
call r3mm(avec,s,s)
! rotate the vectors at each point
lmmax=lmmaxinr
do ir=1,nrmt(is),lrstp
  if (ir.gt.nrmtinr(is)) lmmax=lmmaxvr
  do lm=1,lmmax
    v(:)=srvfmt(lm,ir,:)
    call r3mv(s,v,v)
    srvfmt(lm,ir,:)=v(:)
  end do
end do
return
end subroutine

