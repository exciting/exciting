
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: z2mctm
! !INTERFACE:
subroutine z2mctm(a,b,c)
! !INPUT/OUTPUT PARAMETERS:
!   a   : input matrix 1 (in,complex(2,2))
!   b   : input matrix 2 (in,complex(2,2))
!   c   : output matrix (out,complex(2,2))
! !DESCRIPTION:
!   Multiplies the conjugate transpose of one complex $2\times 2$ matrix with
!   another. Note that the output matrix cannot be one of the input matrices.
!
! !REVISION HISTORY:
!   Created October 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: a(2,2)
complex(8), intent(in) :: b(2,2)
complex(8), intent(out) :: c(2,2)
c(1,1)=conjg(a(1,1))*b(1,1)+conjg(a(2,1))*b(2,1)
c(1,2)=conjg(a(1,1))*b(1,2)+conjg(a(2,1))*b(2,2)
c(2,1)=conjg(a(1,2))*b(1,1)+conjg(a(2,2))*b(2,1)
c(2,2)=conjg(a(1,2))*b(1,2)+conjg(a(2,2))*b(2,2)
return
end subroutine
!EOC


