
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dfqoschd
  implicit none
contains

  subroutine dfqoschd(pou,puo,you,yuo)
    use modmain
    use modxs
    implicit none
    ! arguments
    complex(8), intent(in) :: pou(3),puo(3)
    complex(8), intent(out) :: you,yuo
    you=pou(optcomp(1,1))*conjg(pou(optcomp(2,1)))
    yuo=puo(optcomp(1,1))*conjg(puo(optcomp(2,1)))
  end subroutine dfqoschd

end module m_dfqoschd
