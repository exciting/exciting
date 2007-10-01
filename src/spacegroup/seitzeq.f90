
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

logical function seitzeq(eps,sr1,st1,sr2,st2)
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(in) :: sr1(3,3)
real(8), intent(in) :: st1(3)
real(8), intent(in) :: sr2(3,3)
real(8), intent(in) :: st2(3)
! local variables
integer j
real(8) v1(3),v2(3)
seitzeq=.false.
do j=1,3
  v1(:)=sr1(:,j)+st1(:)
  v2(:)=sr2(:,j)+st2(:)
  if ((abs(v1(1)-v2(1)).gt.eps).or. &
      (abs(v1(2)-v2(2)).gt.eps).or. &
      (abs(v1(3)-v2(3)).gt.eps)) return
end do
seitzeq=.true.
return
end function

