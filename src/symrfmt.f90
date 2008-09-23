
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrfmt
! !INTERFACE:
subroutine symrfmt(lrstp,is,rot,rfmt,srfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   is    : species number (in,integer)
!   rot   : rotation matrix (in,real(3,3))
!   rfmt  : input muffin-tin function (in,real(lmmaxvr,nrmtmax))
!   srfmt : output muffin-tin function (out,real(lmmaxvr,nrmtmax))
! !DESCRIPTION:
!   Applies a symmetry operation (in the form of a rotation matrix) to a real
!   muffin-tin function. The input function can also be the output function. See
!   the routines {\tt rtozflm} and {\tt rotzflm}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
integer, intent(in) :: is
real(8), intent(in) :: rot(3,3)
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax)
real(8), intent(out) :: srfmt(lmmaxvr,nrmtmax)
! local variables
integer ir,irc,nri,nro,iro
! allocatable arrays
complex(8), allocatable :: zfmt1(:,:),zfmt2(:,:)
allocate(zfmt1(lmmaxvr,nrmtmax),zfmt2(lmmaxvr,nrmtmax))
! convert real function to complex spherical harmonic expansion
nri=0
irc=0
do ir=1,nrmt(is),lrstp
  irc=irc+1
  if (ir.le.nrmtinr(is)) then
    call rtozflm(lmaxinr,rfmt(:,ir),zfmt1(:,irc))
    srfmt(lmmaxinr+1:,ir)=0.d0
    nri=irc
  else
    call rtozflm(lmaxvr,rfmt(:,ir),zfmt1(:,irc))
  end if
end do
! first point in the outer point of the muffin-tin
iro=nri+1
! number of points in the outer part
nro=irc-nri
! rotate the complex function
call rotzflm(rot,lmaxinr,nri,lmmaxvr,zfmt1,zfmt2)
call rotzflm(rot,lmaxvr,nro,lmmaxvr,zfmt1(:,iro),zfmt2(:,iro))
! convert complex function to real spherical harmonic expansion
irc=0
do ir=1,nrmt(is),lrstp
  irc=irc+1
  if (ir.le.nrmtinr(is)) then
    call ztorflm(lmaxinr,zfmt2(:,irc),srfmt(:,ir))
  else
    call ztorflm(lmaxvr,zfmt2(:,irc),srfmt(:,ir))
  end if
end do
deallocate(zfmt1,zfmt2)
return
end subroutine
!EOC

