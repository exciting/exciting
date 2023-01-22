!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writestate
! !INTERFACE:
!
!
Subroutine writestate
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Writes the charge density, potentials and other relevant variables to the
!   file {\tt STATE.OUT}. Note to developers: changes to the way the variables
!   are written should be mirrored in {\tt readstate}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is
      Open (50, File='STATE'//trim(filext), Action='WRITE', Form='UNFOR&
     &MATTED')
      Write (50) version, githash
      Write (50) associated (input%groundstate%spin)
      Write (50) nspecies
      Write (50) lmmaxvr
      Write (50) nrmtmax
      Do is = 1, nspecies
         Write (50) natoms (is)
         Write (50) nrmt (is)
         Write (50) spr (1:nrmt(is), is)
      End Do
      Write (50) ngrid
      Write (50) ngvec
      Write (50) ndmag
      Write (50) nspinor
      Write (50) ldapu
      Write (50) lmmaxlu
! write the density
      Write (50) rhomt, rhoir
! write the Coulomb potential
      Write (50) vclmt, vclir
! write the exchange-correlation potential
      Write (50) vxcmt, vxcir
! write the effective potential
      Write (50) veffmt, veffir, veffig
! write the magnetisation and effective magnetic fields
      If (associated(input%groundstate%spin)) Then
         Write (50) magmt, magir
         Write (50) bxcmt, bxcir
      End If
! write the LDA+U potential matrix elements
      If (ldapu .Ne. 0) Then
         Write (50) vmatlu
      End If
      Close (50)
      Return
End Subroutine
!EOC
