!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine writechg (fnum)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: fnum
! local variables
      Integer :: is, ia, ias
! output charges
      Write (fnum,*)
      Write (fnum, '("Charges :")')
      Write (fnum, '(" core", T30, ": ", G18.10)') chgcr
      Write (fnum, '(" core leakage", T30, ": ", G18.10)') chgcrlk
      Write (fnum, '(" valence", T30, ": ", G18.10)') chgval
      Write (fnum, '(" interstitial", T30, ": ", G18.10)') chgir
      Write (fnum, '(" muffin-tins")')
      Do is = 1, nspecies
         Write (fnum, '("  species : ", I4, " (", A, ")")') is, trim &
        & (input%structure%speciesarray(is)%species%chemicalSymbol)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Write (fnum, '("   atom ", I4, T30, ": ", G18.10)') ia, &
           & chgmt (ias)
         End Do
      End Do
      Write (fnum, '(" total in muffin-tins", T30, ": ", G18.10)') &
     & chgmttot
      If (input%groundstate%chgexs .Ne. 0.d0) Then
         Write (fnum, '(" excess", T30, ": ", G18.10)') &
        & input%groundstate%chgexs
      End If
      Write (fnum, '(" total charge", T30, ": ", G18.10)') chgcalc
! output moments
      If (associated(input%groundstate%spin)) Then
         Write (fnum,*)
         Write (fnum, '("Moments :")')
         Write (fnum, '(" interstitial", T30, ": ", 3G18.10)') momir &
        & (1:ndmag)
         Write (fnum, '(" muffin-tins")')
         Do is = 1, nspecies
            Write (fnum, '("  species : ", I4, " (", A, ")")') is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol)
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Write (fnum, '("	atom ", I4, T30, ": ", 3G18.10)') ia, &
              & mommt (1:ndmag, ias)
            End Do
         End Do
         Write (fnum, '(" total in muffin-tins", T30, ": ", 3G18.10)') &
        & mommttot (1:ndmag)
         Write (fnum, '(" total moment", T30, ": ", 3G18.10)') momtot &
        & (1:ndmag)
      End If
!call flushifc(fnum)
      Return
End Subroutine
