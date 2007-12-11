
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dfqoscbo
  implicit none
contains

  subroutine dfqoscbo(n,xou,xuo,you,yuo)
    use m_tdzoutpr2
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: xou(:),xuo(:)
    complex(8), intent(out) :: you(:,:),yuo(:,:)
    ! local variables
    complex(8) :: zt
    zt=(1.d0,0.d0)
    call tdzoutpr2(n,n,zt,xou,xou,you)
    call tdzoutpr2(n,n,zt,xuo,xuo,yuo)
  end subroutine dfqoscbo

end module m_dfqoscbo
