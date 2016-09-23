!
!
!
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematqk1 (iq, ik)
      Use modmain
      Use modinput
      Use modxs
      Use modmpi
      Use m_putemat
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik
  ! set band combinations
      If ( .Not. (task .Eq. 430)) Then
         Call ematbdlims (2*input%xs%emattype, nst1, istl1, istu1, &
        & nst2, istl2, istu2)
         If (allocated(xiou)) deallocate (xiou)
         If (allocated(xiuo)) deallocate (xiuo)
         Allocate (xiou(nst1, nst2, ngq(iq)))
         Call ematqk (iq, ik)
      End If
      If (input%xs%emattype .Eq. 0) Then
     ! all band combinations
         nst3 = nstsv
         nst4 = nstsv
         If ( .Not. ((task .Ge. 400) .And. (task .Le. 499))) Call &
        & putemat (iq, ik, .True., trim(fnemat), istl1, istu1, istl2, &
        & istu2, xiou)
      Else
     ! o-u/u-o or o-o/u-u band combinations
         If ( .Not. (task .Eq. 430)) Then
            Allocate (xiuo(nst1, nst2, ngq(iq)))
            xiuo (:, :, :) = xiou (:, :, :)
         End If
         Call ematbdlims (2*input%xs%emattype-1, nst1, istl1, istu1, &
        & nst2, istl2, istu2)
         istl3 = istl2
         istu3 = istu2
         istl4 = istl1
         istu4 = istu1
         nst3 = nst2
         nst4 = nst1
         If (allocated(xiou)) deallocate (xiou)
         Allocate (xiou(nst1, nst2, ngq(iq)))
         Call ematqk (iq, ik)
         If ( .Not. tscreen) Call putemat (iq, ik, .True., &
        & trim(fnemat), istl1, istu1, istl2, istu2, xiou, istl3, istu3, &
        & istl4, istu4, xiuo)
      End If
End Subroutine ematqk1
