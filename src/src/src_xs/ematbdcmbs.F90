!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematbdcmbs (etyp)
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: etyp
      Select Case (etyp)
      Case (0)
     ! all combinations
         Call ematbdlims (0, nst1, istl1, istu1, nst2, istl2, istu2)
         nst3 = 0
         nst4 = 0
      Case (1)
     ! o-u combinations
         Call ematbdlims (1, nst1, istl1, istu1, nst2, istl2, istu2)
     ! u-o combinations
         Call ematbdlims (2, nst3, istl3, istu3, nst4, istl4, istu4)
      Case (2)
     ! o-o combinations
         Call ematbdlims (3, nst1, istl1, istu1, nst2, istl2, istu2)
     ! u-u combinations
         Call ematbdlims (4, nst3, istl3, istu3, nst4, istl4, istu4)
      End Select
End Subroutine ematbdcmbs
