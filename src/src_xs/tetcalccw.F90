!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine tetcalccw
      Use modmain
      Use modinput
      Use modxs
      Use modmpi
      Use modtetra
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'tetcalccw'
      Integer :: iq
      Logical :: tet
      Call init0
  ! initialise universal variables
      tet = input%xs%tetra%tetradf
      input%xs%tetra%tetradf = .True.
      Call init1
  ! save Gamma-point variables
      Call xssave0
  ! initialize q-point set
      Call init2
  ! read Fermi energy
      Call readfermi
  ! w-point interval for process
      If (tscreen) Then
         nwdf = 1
         Call genparidxran ('q', nqpt)
      Else
         Call genparidxran ('w', nwdf)
      End If
  ! loop over q-points
      Do iq = qpari, qparf
     ! call for q-point
         Call tetcalccwq (iq)
         Write (unitout, '(a, i8)') 'Info(' // thisnam // '): weights f&
        &or tetrahedron method finished for q - point:', iq
      End Do
  ! synchronize
      Call barrier
      If ((procs .Gt. 1) .And. (rank .Eq. 0) .And. ( .Not. tscreen)) &
     & Call tetgather
      Call barrier
      input%xs%tetra%tetradf = tet
      Write (unitout, '(a)') "Info(" // trim (thisnam) // "): weights f&
     &or tetrahedron method finished"
      Call genfilname (setfilext=.True.)
End Subroutine tetcalccw
