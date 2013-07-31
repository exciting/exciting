!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine writechg (fnum,verbosity)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: fnum
      integer, intent(in) :: verbosity
! local variables
      Integer :: is, ia, ias
      character*(77) :: string

!_______________
! output charges

      if ( verbosity < 1 ) then
          Write (fnum, '(" Electron charge: Core leakage", T45, ": ", F18.8)') chgcrlk
          return
      end if 

      Write (fnum,*)
      Write (fnum, '(" Electron charges :")')
      Write (fnum, '("     core", T45, ": ", F18.8)') chgcr
      Write (fnum, '("     core leakage", T45, ": ", F18.8)') chgcrlk
      Write (fnum, '("     valence", T45, ": ", F18.8)') chgval
      Write (fnum, '("     interstitial", T45, ": ", F18.8)') chgir
      Write (fnum, '("     charge in muffin-tin spheres :")')
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            write(fnum,'("                  atom ",I5,4x,A2,T45,": ",F18.8)') &
           &  ia, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
           &  chgmt (ias)
         End Do
      End Do
      Write (fnum, '("     total charge in muffin-tins", T45, ": ", F18.8)') &
     & chgmttot
      If (input%groundstate%chgexs .Ne. 0.d0) Then
         Write (fnum, '("     excess charge", T45, ": ", F18.8)') &
        & input%groundstate%chgexs
      End If
      Write (fnum, '("     total charge", T45, ": ", F18.8)') chgcalc

!_______________
! output moments

      If (associated(input%groundstate%spin)) Then
         Write (fnum,*)
         Write (fnum, '(" Moments :")')
         Write (fnum, '("     interstitial", T28, ": ", 3F16.8)') momir &
        & (1:ndmag)
         Write (fnum, '("     muffin-tins")')
         Do is = 1, nspecies
            Write (fnum, '("     species : ", I4, " (", A, ")")') is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol)
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Write (fnum, '("     atom ", I4, T28, ": ", 3F16.8)') ia, &
              & mommt (1:ndmag, ias)
            End Do
         End Do
         Write (fnum, '("     total in muffin-tins", T28, ": ", 3F16.8)') &
        & mommttot (1:ndmag)
         Write (fnum, '("     total moment", T28, ": ", 3F16.8)') momtot &
        & (1:ndmag)
      End If

      call flushifc(fnum)

      Return
End Subroutine
