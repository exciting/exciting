
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gndstateq
  implicit none
contains

  subroutine gndstateq(voff,filxt)
    use modmain
    implicit none
    ! arguments
    real(8), intent(in) :: voff(:)
    character(*), intent(in) :: filxt
    ! local varialbes
    character(256) :: filext_save
    real(8) :: vkloff_save(3)
    integer :: task_save, maxscl_save

    ! save original values
    filext_save=trim(filext)
    vkloff_save=vkloff
    task_save=task
    maxscl_save=maxscl

    ! one iteration, new offset, special file extension
    filext=trim(filxt)
    vkloff=voff
    task=1
    maxscl=1
    
    ! call with the above parameters changed
    call gndstate
    
    ! restore original parameters
    filext=trim(filext_save)
    vkloff=vkloff_save
    task=task_save
    maxscl=maxscl_save
    
  end subroutine gndstateq

end module m_gndstateq
