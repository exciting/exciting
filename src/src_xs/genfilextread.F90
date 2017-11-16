! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genfilextread
! !INTERFACE:
subroutine genfilextread(task)
! !USES:
  use modinput
  use m_genfilname
! !INPUT/OUTPUT PARAMETERS:
! !IN:
! integer(4) :: task ! Excitings faboulus task number
! 
! !DESCRIPTION:
!   This routine sets the filename extension for the
!   eigenvalue {\tt EVALSV} and eigenvector {\tt EVECSV} files
!   corresponding to zero momentum transfer depending on the calling task.
!   Used in {\tt findocclims}.
!
! !REVISION HISTORY:
!   Added to documentaion scheme. 2016 (Aurich)
!   Changed so that only task 'screen' uses '_SCR.OUT' (Aurich)
!
!EOP
!BOC
  implicit none

  integer(4), intent(in) :: task

  select case(task)
    ! Task='screen'
    case(430)
      call genfilname(iqmt=1, scrtype='', setfilext=.true.)
    case default
      call genfilname(iqmt=1, setfilext=.true.)
  end select

end subroutine genfilextread
!EOC
