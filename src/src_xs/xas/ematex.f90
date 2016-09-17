! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ematex
! !INTERFACE:
Subroutine ematex(iq, ik)
! !USES:
      Use modmain
      Use modinput
      Use modxs
      Use modmpi
      Use m_putemat
      Use modxas
! !INPUT/OUTPUT PARAMETERS:
!   iq       : q-point position (in,integer)
!   ik       : k-point position (in,integer)
! !DESCRIPTION:
!   Calculates Planewave Matrix elements for the exchange term of the BSE Hamiltonian.
!
! !REVISION HISTORY:
!  Based on the subroutine ematqk1.F90
!  Created November 2015 (Christian Vorwerk)
!EOP
!BOC      

      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik
      Integer :: n1, n2
  ! set band combinations
      If ( .Not. (task .Eq. 430)) Then
         Call ematbdlims (2, nst1, istl1, istu1, &
        & nst2, istl2, istu2)
         If (allocated(xiou)) deallocate (xiou)
         If (allocated(xiuo)) deallocate (xiuo)
         !Allocate (xiuo(nst1, nst2, ngq(iq)))
         !Call genxiuo (iq, ik, xiuo)
      End If
      If (input%xs%emattype .Eq. 0) Then
     ! all band combinations (only for debugging!)
         nst3 = nstsv
         nst4 = nstsv
         If ( .Not. ((task .Ge. 400) .And. (task .Le. 499))) Call &
        & putemat (iq, ik, .True., trim(fnemat), istl1, istu1, istl2, &
        & istu2, xiuo)
      Else
     ! o-u/u-o or o-o/u-u band combinations
         Call ematbdlims (1, nst1, istl1, istu1, &
        & nst2, istl2, istu2)
         istl3 = istl2
         istu3 = istu2
         istl4 = istl1
         istu4 = istu1
         nst3 = nst2
         nst4 = nst1
         If (allocated(xiou)) deallocate (xiou)
         Allocate (xiou(nst1, nst2, ngq(iq)))
         ! Call to Subroutine genxiou
         Call genxiou (iq, ik, xiou)
         If ( .Not. tscreen) Call putemat (iq, ik, .True., &
        & trim(fnemat), istl1, istu1, istl2, istu2, xiuo, istl3, istu3, &
        & istl4, istu4, xiou)
      End If
End Subroutine ematex
!EOC
