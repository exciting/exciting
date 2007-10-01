
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: radmesh
! !INTERFACE:
subroutine radmesh(nr,irf,rf,rmin,r)
! !INPUT/OUTPUT PARAMETERS:
!   nr   : number of mesh points (in,integer)
!   irf  : fixed point number (in,integer)
!   rf   : fixed point radius (in,real)
!   rmin : minimum radius (in,real)
!   r    : radial mesh (out,real(n))
! !DESCRIPTION:
!   Generates a logarithmic radial mesh containing $n_r$ points. The mesh
!   contains a fixed point $i_{\rm f}$ corresponding to radius $r_{\rm f}$.
!   Mesh radii are given by
!   $$ r(i)=r_{\rm min}\exp\left(\frac{i-1}{i_{\rm f}-1}t\right), $$
!   where $t=\log(r_{\rm f}/r_{\rm min})$.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr
integer, intent(in) :: irf
real(8), intent(in) :: rf
real(8), intent(in) :: rmin
real(8), intent(out) :: r(nr)
! local variables
integer ir
real(8) x,t1,t2
if (nr.lt.2) then
  write(*,*)
  write(*,'("Error(radmesh): nr < 2 : ",I8)') nr
  write(*,*)
  stop
end if
if (rmin.lt.1.d-20) then
  write(*,*)
  write(*,'("Error(radmesh): rmin too small : ",G18.10)') rmin
  write(*,*)
  stop
end if
if ((irf.lt.2).or.(irf.gt.nr)) then
  write(*,*)
  write(*,'("Error(radmesh): invalid irf : ",I8)') irf
  write(*,*)
  stop
end if
t1=log(rf/rmin)
t2=1.d0/dble(irf-1)
do ir=1,nr
  x=dble(ir-1)*t2
  r(ir)=rmin*exp(x*t1)
end do
return
end subroutine
!EOC

