
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dfqoscbo
  implicit none
contains

  subroutine dfqoscbo(n,xou,xuo,you,yuo)
    use m_xszoutpr2
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: xou(:),xuo(:)
    complex(8), intent(out) :: you(:,:),yuo(:,:)
    ! local variables
    complex(8) :: zt

    zt=(1.d0,0.d0)
    call xszoutpr2(n,n,zt,xou,xou,you)
    call xszoutpr2(n,n,zt,xuo,xuo,yuo)

  end subroutine dfqoscbo

end module m_dfqoscbo
