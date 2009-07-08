

! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine seitzmul(eps, sr1, st1, sr2, st2, sr3, st3)
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(in) :: sr1(3, 3)
real(8), intent(in) :: st1(3)
real(8), intent(in) :: sr2(3, 3)
real(8), intent(in) :: st2(3)
real(8), intent(out) :: sr3(3, 3)
real(8), intent(out) :: st3(3)
! local variables
integer::id(3)
call r3mv(sr1, st2, st3)
st3(:)=st3(:)+st1(:)
call r3frac(eps, st3, id)
call r3mm(sr1, sr2, sr3)
return
end subroutine
