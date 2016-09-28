!
!
!
! Copyright (C) 2015 Christian Vorwerk and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
!
! !ROUTINE: ematdirect
!
! !INTERFACE:
Subroutine ematdirect (iq, ik)
! !INPUT/OUTPUT PARAMETERS:
!   iq       : q-point position (in,integer)
!   ik       : k-point position (in,integer)
! !DESCRIPTION:
!
! Calculates the plane wave matrix elements $$ M_{cc' \bf k}(\bf q + \bf G)$$ between two core states and  $$ M_{uu' \bf k}(\bf q + \bf G)$$ between two conduction states.
! The matrix elements are used for the construction of the direct terms in the BSE Hamiltonian (see xas.f90).
! 
! !USES:
      Use modmain
      Use modinput
      Use modxs
      Use modmpi
      Use m_putemat
      Use modxas
      Implicit None
! !REVISION HISTORY:
!
! Created June 2015 by C. Vorwerk
!
!EOP
!BOC
  ! arguments
      Integer, Intent (In) :: iq, ik
  ! set band combinations
      If ( .Not. (task .Eq. 430)) Then
         Call ematbdlims (4, nst1, istl1, istu1, &
        & nst2, istl2, istu2)
         If (allocated(xiou)) deallocate (xiou)
         If (allocated(xiuo)) deallocate (xiuo)
         Allocate (xiou(nst1, nst2, ngq(iq)))
         ! Calculation of Matrix elements between two conduction states
         Call ematqk (iq, ik)
      End If
     ! o-u/u-o or o-o/u-u band combinations
         If ( .Not. (task .Eq. 430)) Then
			If (allocated(xiuu)) deallocate (xiuu)
            Allocate (xiuu(nst1, nst2, ngq(iq)))
            xiuu (:, :, :) = xiou (:, :, :)
         End If
         Call ematbdlims (3, nst1, istl1, istu1, &
        & nst2, istl2, istu2)
         istl3 = istl2
         istu3 = istu2
         istl4 = istl1
         istu4 = istu1
         nst3 = nst2
         nst4 = nst1
         If (allocated(xioo)) deallocate (xioo)
         Allocate (xioo(nst1, nst2, ngq(iq)))
         ! Calculation of matrix elements between two core states
         Call genxioo (iq, ik, xioo)
         If ( .Not. tscreen) Call putemat (iq, ik, .True., &
        & trim(fnemat), istl1, istu1, istl2, istu2, xiuu, istl3, istu3, &
        & istl4, istu4, xioo)
End Subroutine ematdirect
! EOC
