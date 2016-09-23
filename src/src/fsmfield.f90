!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: fsmfield
! !INTERFACE:
!
!
Subroutine fsmfield
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Updates the effective magnetic field, ${\bf B}_{\rm FSM}$, required for
!   fixing the spin moment to a given value, $\boldsymbol{\mu}_{\rm FSM}$. This
!   is done by adding a vector to the field which is proportional to the
!   difference between the moment calculated in the $i$th self-consistent loop
!   and the required moment:
!   $$ {\bf B}_{\rm FSM}^{i+1}={\bf B}_{\rm FSM}^i+\lambda\left(
!    \boldsymbol{\mu}^i-\boldsymbol{\mu}_{\rm FSM}\right), $$
!   where $\lambda$ is a scaling factor.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, idm
      Real (8) :: v (3), t1
      If (( .Not. associated(input%groundstate%spin)) .Or. &
     & (input%groundstate%spin%fixspinnumber .Eq. 0)) Return
      t1 = 1.d0 / y00
! determine the global effective field
      If ((input%groundstate%spin%fixspinnumber .Eq. 1) .Or. &
     & (input%groundstate%spin%fixspinnumber .Eq. 3)) Then
         If (ncmag) Then
            v (:) = input%groundstate%spin%momfix(:)
         Else
            v (1) = input%groundstate%spin%momfix(3)
         End If
         Do idm = 1, ndmag
            bfsmc (idm) = bfsmc (idm) + input%groundstate%spin%taufsm * &
           & (momtot(idm)-v(idm))
         End Do
         Do idm = 1, ndmag
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  Do ir = 1, nrmt (is)
                     bxcmt (1, ir, ias, idm) = bxcmt (1, ir, ias, idm) &
                    & + t1 * bfsmc (idm)
                  End Do
               End Do
            End Do
            Do ir = 1, ngrtot
               bxcir (ir, idm) = bxcir (ir, idm) + bfsmc (idm)
            End Do
         End Do
      End If
      If ((input%groundstate%spin%fixspinnumber .Eq. 2) .Or. &
     & (input%groundstate%spin%fixspinnumber .Eq. 3)) Then
! determine the muffin-tin fields for fixed local moments
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               If (ncmag) Then
                  v (:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%mommtfix(:)
               Else
                  v (1) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%mommtfix(3)
               End If
               Do idm = 1, ndmag
                  bfsmcmt (idm, ia, is) = bfsmcmt (idm, ia, is) + &
                 & input%groundstate%spin%taufsm * (mommt(idm, &
                 & ias)-v(idm))
                  Do ir = 1, nrmt (is)
                     bxcmt (1, ir, ias, idm) = bxcmt (1, ir, ias, idm) &
                    & + t1 * bfsmcmt (idm, ia, is)
                  End Do
               End Do
            End Do
         End Do
      End If
      Return
End Subroutine
!EOC
