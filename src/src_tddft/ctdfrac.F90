
! Copyright (C) 2006 S. Sagmeister.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ctdfrac
  implicit none
contains

  subroutine ctdfrac(n,a,b,f)
    ! Straight foreward implemetation of a continued fraction with depth n
    ! without any checking of convergence.
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: a(n), b(0:n)
    complex(8), intent(out) :: f  
    ! local variables
    integer :: i,j
    f=b(n)
    do j=n,1,-1
      f=b(j-1)+a(j)/f
    end do
  end subroutine

end module

