!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writedveff (iq, is, ia, ip, dveffmt, dveffir)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iq
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ip
      Complex (8), Intent (In) :: dveffmt (lmmaxvr, nrcmtmax, natmtot0)
      Complex (8), Intent (In) :: dveffir (ngrtot0)
! local variables
      Integer :: js
      Character (256) :: fext
      Call phfext (iq, is, ia, ip, fext)
      Open (50, File='DVEFF'//trim(fext), Action='WRITE', Form='UNFORMA&
     &TTED')
      Write (50) version
      Write (50) nspecies
      Write (50) lmmaxvr
      Do js = 1, nspecies
         Write (50) natoms0 (js)
         Write (50) nrcmt (js)
      End Do
      Write (50) ngrid0
      Write (50) dveffmt, dveffir
      Close (50)
      Return
End Subroutine
