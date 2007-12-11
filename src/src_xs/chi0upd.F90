
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_chi0upd
  implicit none
contains

  subroutine chi0upd(n,wou,wuo,hou,huo,chi0)
    use modmain
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: wou,wuo,hou(n,n),huo(n,n)
    complex(8), intent(inout) :: chi0(n,n)
    ! local variables
    integer :: i,m
    real(8), allocatable :: bhou(:,:), chou(:,:), bhuo(:,:), chuo(:,:)
    real(8), allocatable :: bchi0(:,:),cchi0(:,:)

    do m=1,n
       call zaxpy(n,wou,hou(1,m),1,chi0(1,m),1)
       call zaxpy(n,wuo,huo(1,m),1,chi0(1,m),1)
    end do

  end subroutine chi0upd

end module m_chi0upd
