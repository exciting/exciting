!
!
!
! Copyright (C) 2014 S. Kontur and C. Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genkcpts
!
Subroutine bsedgridinit ()
! !USES:
     Use modmain
     Use modxs
     Use modinput
     Use m_gndstateq
     Use modmpi
     Use m_genfilname
     Use m_getunit
     Implicit None
     Character (77) :: string
! !DESCRIPTION:
!   Initializes for one BSE calculation within a double grid loop.
!
! !REVISION HISTORY:
!   Created January 2014, S. Kontur
!EOP
!BOC
 
     ! add subgrid offset to the one requested in input file for XS
     input%xs%vkloff = vksubl(:, iksubpt) + vkloff_xs_b(:)
     input%groundstate%vkloff = input%xs%vkloff
     input%xs%screening%vkloff = input%xs%vkloff
     input%xs%BSE%vkloff = input%xs%vkloff
     input%xs%screening%reducek = .false.
     input%xs%BSE%reducek = .false.
     if (allocated(evalsv0)) deallocate (evalsv0)
!    call xsinit
!    call init0
!    call init1

     Call genfilname (nodotpar=.True., basename='INFOXS', procs=procs, rank=rank, filnam=xsfileout)
     Call getunit (unitout)
     open (unit=unitout, file=xsfileout, status="unknown", action="write", position="append")
     write (string,'("Sub grid point ",i3," with offset ",3f7.3)') iksubpt, vksubl(:, iksubpt)+vkloff_xs_b(:)
     call printline (unitout,"=")
     call printtext (unitout,"=",string)
     call printline (unitout,"=")
     close (unitout)

     call gndstateq (input%xs%vkloff, '.OUT')

End Subroutine bsedgridinit
!EOC
