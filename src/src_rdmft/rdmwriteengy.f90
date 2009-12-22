!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmwriteengy (fnum)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: fnum
      Write (fnum,*)
      Write (fnum, '("Energies :")')
      Write (fnum, '(" electronic kinetic", T30, ": ", G18.10)') engykn
      Write (fnum, '(" core electron kinetic", T30, ": ", G18.10)') &
     & engykncr
      Write (fnum, '(" Coulomb", T30, ": ", G18.10)') engyvcl
      Write (fnum, '(" Madelung", T30, ": ", G18.10)') engymad
      Write (fnum, '(" exchange-correlation", T30, ": ", G18.10)') &
     & engyx
      If (input%groundstate%RDMFT%rdmtemp .Gt. 0.d0) Then
         Write (*, '(" entropy", T30, ": ", G18.10)') rdmentrpy
      End If
      Write (fnum, '(" total energy", T30, ": ", G18.10)') engytot
      Call flushifc (fnum)
      Return
End Subroutine
