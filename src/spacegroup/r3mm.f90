
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3mm
! !INTERFACE:
subroutine r3mm(a,b,c)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix 1 (in,real(3,3))
!   b : input matrix 2 (in,real(3,3))
!   c : output matrix (out,real(3,3))
! !DESCRIPTION:
!   Multiplies two real $3\times 3$ matrices. The output matrix can also be an
!   input matrix.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: a(3,3)
real(8), intent(in) :: b(3,3)
real(8), intent(out) :: c(3,3)
! local variables
real(8) d(3,3)
d(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
d(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)+a(2,3)*b(3,1)
d(3,1)=a(3,1)*b(1,1)+a(3,2)*b(2,1)+a(3,3)*b(3,1)
d(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
d(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
d(3,2)=a(3,1)*b(1,2)+a(3,2)*b(2,2)+a(3,3)*b(3,2)
d(1,3)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
d(2,3)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
d(3,3)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
! copy to output matrix
c(:,:)=d(:,:)
return
end subroutine
!EOC
