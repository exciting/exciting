
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine grad2rfmt(lmax,nr,r,ld,rfmt,g2rfmt)
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld,nr)
real(8), intent(out) :: g2rfmt(ld,nr)
! local variables
integer l,m,lm,ir
real(8) t1,f(nr)
! automatic arrays
real(8) ri(nr),ri2(nr),cf(4,nr)
! tabulate 1/r and 1/r^2
do ir=1,nr
  ri(ir)=1.d0/r(ir)
  ri2(ir)=ri(ir)**2
end do
lm=0
do l=0,lmax
  t1=-dble(l*(l+1))
  do m=-l,l
    lm=lm+1
! use a cubic spline to compute radial derivatives
   do ir=1,nr
     f(ir)=r(ir)*rfmt(lm,ir)
   end do
!    call spline4(nr,r,ld,rfmt(lm,1),cf)
    call spline4(nr,r,1,f,cf)
! apply Laplacian
    do ir=1,nr
!      g2rfmt(lm,ir)=2.d0*(ri(ir)*cf(1,ir)+cf(2,ir))+ri2(ir)*t1*rfmt(lm,ir)
      g2rfmt(lm,ir)=2.d0*cf(2,ir)*ri(ir)+ri2(ir)*t1*rfmt(lm,ir)
    end do
  end do
end do


 return
end subroutine

