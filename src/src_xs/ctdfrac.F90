
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ctdfrac
  implicit none
contains

!BOP
! !ROUTINE: ctdfrac
! !INTERFACE:
  subroutine ctdfrac(n,a,b,f)
! !INPUT/OUTPUT PARAMETERS:
!   n     : depth of continued fraction (in,integer)
!   a     : a-coefficients (in,complex(n))
!   b     : b-coefficients (in,complex(0:n))
!   f     : continued fraction result
! !DESCRIPTION:
!   Straight foreward evaluation of a continued fraction with depth $n$
!   $$ b_0+\cfrac{a_1}{b_1+\cfrac{a_2}
!   {b_2+\cfrac{a_3}{\cdots+\cfrac{a_n}{b_n}}}}  $$
!   without any checking of convergence.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: a(n), b(0:n)
    complex(8), intent(out) :: f  
    ! local variables
    integer :: j
    f=b(n)
    do j=n,1,-1
      f=b(j-1)+a(j)/f
    end do
  end subroutine

end module
!EOC
