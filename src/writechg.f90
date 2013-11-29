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
      Integer :: is, ia, ias, i
      character*(77) :: string

!_______________
! output charges

      if ( verbosity < 1 ) then
          Write (fnum, '(" Electron charge: Core leakage", T45, ": ", F18.8)') chgcrlk
          return
      end if 

      i = 0
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
            i = i+1
            write(fnum,'("                  atom ",I5,4x,A2,T45,": ",F18.8)') &
           &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), chgmt (ias)
         End Do
      End Do
      Write (fnum, '("     total charge in muffin-tins", T45, ": ", F18.8)') chgmttot
      If (input%groundstate%chgexs .Ne. 0.d0) Then
         Write (fnum, '("     excess charge", T45, ": ", F18.8)') input%groundstate%chgexs
      End If
      Write (fnum, '("     total charge", T45, ": ", F18.8)') chgcalc

!_______________
! output moments

      If (associated(input%groundstate%spin)) Then
         i = 0
         Write (fnum,*)
         Write (fnum, '(" Moments :")')
         if (ndmag==1) then
            Write (fnum, '("     interstitial", T45, ": ", F18.8)') momir (1:ndmag)
            Write (fnum, '("     moment in muffin-tin spheres :")')
            do is = 1, nspecies
               do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  i = i+1
                  write(fnum,'("                  atom ",I5,4x,A2,T45,": ",F18.8)') &
                 &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), mommt(1:ndmag, ias)
               end do
            end do
            Write (fnum, '("     total moment in muffin-tins", T45, ": ", F18.8)') mommttot (1:ndmag)
            Write (fnum, '("     total moment", T45, ": ", F18.8)') momtot (1:ndmag)
         else
            Write (fnum, '("     interstitial", T35, ":", 3F15.8)') momir (1:ndmag)
            Write (fnum, '("     moment in muffin-tin spheres :")')
            do is = 1, nspecies
               do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  i = i+1
                  write(fnum,'("         atom ",I5,4x,A2,T35,":",3F15.8)') &
                 &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), mommt(1:ndmag, ias)
               end do
            end do
            Write (fnum, '("     total moment in muffin-tins", T35, ":", 3F15.8)') mommttot (1:ndmag)
            Write (fnum, '("     total moment", T35, ":", 3F15.8)') momtot (1:ndmag)
         end if

      End If

      call flushifc(fnum)

      Return
End Subroutine
