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
Subroutine genksubpts ()
! !USES:
     Use modmain
     Use modxs
     Use modinput
     Use m_getunit
     Use m_genfilname
     Use modmpi
! !DESCRIPTION:
!   Sets up a double k-grid for BSE calculations.
!
! !REVISION HISTORY:
!   Created January 2014, S. Kontur
!EOP
!BOC
     Implicit None
! local variables
     Integer :: un, ngrids
     Real (8) :: boxl (3, 4)
     Character (77) :: string

     ngrids = product(input%xs%BSE%ngridksub(1:3))

! allocate the reduced k-point set arrays
     If (allocated(ivksub)) deallocate (ivksub)
     Allocate (ivksub(3, ngrids))
     If (allocated(vksubl)) deallocate (vksubl)
     Allocate (vksubl(3, ngrids))
     If (allocated(vksubc)) deallocate (vksubc)
     Allocate (vksubc(3, ngrids))
     If (allocated(wksubpt)) deallocate (wksubpt)
     Allocate (wksubpt(input%xs%BSE%ngridksub(1)*&
        & input%xs%BSE%ngridksub(2)*input%xs%BSE%ngridksub(3)))
     If (allocated(iksubmap)) deallocate (iksubmap)
     Allocate (iksubmap(0:input%xs%BSE%ngridksub(1)-1, &
        & 0:input%xs%BSE%ngridksub(2)-1, &
        & 0:input%xs%BSE%ngridksub(3)-1))

! setup the default k-point box (without offset)
     boxl (:, 1) = dble (input%xs%BSE%ngridksub(:))
     boxl (:, 2) = boxl (:, 1)
     boxl (:, 3) = boxl (:, 1)
     boxl (:, 4) = boxl (:, 1)
     boxl (1, 2) = boxl (1, 2) + 1.d0
     boxl (2, 3) = boxl (2, 3) + 1.d0
     boxl (3, 4) = boxl (3, 4) + 1.d0
     Call genppts (.True., .False., input%xs%BSE%ngridksub, boxl,  &
  &        nksubpt, iksubmap, ivksub, vksubl, vksubc, wksubpt)
     Call getunit (un)
     Open (un, File='KSUBPOINTS'//trim(filext), Action='WRITE', Form='FORMATTED')
     Write (un, '(I6, " : nksubpt; sub k-point, vksubl, wksubpt below")') nkpt
     Do iksubpt = 1, nksubpt
        Write (un, '(I6, 4G18.10, 2I8)') iksubpt, vksubl (:, iksubpt), wksubpt (iksubpt)
     End Do
     Close (un)

     Call genfilname (nodotpar=.True., basename='INFOXS', procs=procs, rank=rank, filnam=xsfileout)
     Call getunit (unitout)
     open (unit=unitout, file=xsfileout, status="unknown", action="write")
     write (string,'("BSE double grid computation: subgrid generated")')
     call printline (unitout,"=")
     call printtext(unitout,"=",string)
     call printline (unitout,"=")
     close (unitout)

End Subroutine genksubpts
!EOC
