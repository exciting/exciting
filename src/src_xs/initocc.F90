!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine initocc (nbf, naf)
      Use modmain
      Use modinput
      Use modxs
      Implicit None
  ! arguments
      Integer :: nbf, naf
  ! local variables
      Character (*), Parameter :: thisnam = 'initocc'
      Integer :: nvalel
  ! number of valence electrons
      nvalel = Nint (chgval/2.d0)
  ! number of states below Fermi energy
      If (nbf .Eq. 0) nbf = nvalel
      If (nbf .Gt. nvalel) Then
         Write (unitout, '("Warning(", a, "): number of states below Fe&
        &rmi energy too large - adjusting to number of valence states")&
        &') trim (thisnam)
      End If
  ! number of states above Fermi energy
      If (naf .Eq. 0) Then
     ! if "naf" is not specified define it using "nempty"
         naf = input%xs%nempty + 1
      Else
     ! check if number is too large
         If (naf .Gt. (input%xs%nempty+1)) Then
            naf = input%xs%nempty + 1
            Write (unitout, '("Warning(", a, "): number of states above&
           & Fermi energy too large - adjusting using number of empty s&
           &tates")') trim (thisnam)
         End If
      End If
End Subroutine initocc
