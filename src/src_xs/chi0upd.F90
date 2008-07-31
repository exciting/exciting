
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_chi0upd
  implicit none
contains

!BOP
! !ROUTINE: chi0upd
! !INTERFACE:
  subroutine chi0upd(n,wou,wuo,hou,huo,chi0)
! !USES:
    use modmain
! !DESCRIPTION:
!   Updates the Kohn-Sham response function.
!
! !REVISION HISTORY:
!   Created January 2005 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: wou,wuo,hou(n,n),huo(n,n)
    complex(8), intent(inout) :: chi0(n,n)
    ! local variables
    integer :: m
    do m=1,n
       call zaxpy(n,wou,hou(1,m),1,chi0(1,m),1)
       call zaxpy(n,wuo,huo(1,m),1,chi0(1,m),1)
    end do
  end subroutine chi0upd

end module m_chi0upd
!EOC
