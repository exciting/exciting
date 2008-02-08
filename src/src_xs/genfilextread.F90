 
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genfilextread(task)
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: task
  select case(task)
  case(330,340,121)
     call genfilname(iqmt=0,setfilext=.true.)
  case(430)
     call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  case(440)
     call genfilname(dotext='_SCI.OUT',setfilext=.true.)
  end select
end subroutine genfilextread
