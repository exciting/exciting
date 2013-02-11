!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: charge
! !INTERFACE:
!
!
Subroutine charge
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total charges by integrating the
!   density.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir
      Real (8) :: sum, t1
! automatic arrays
      Real (8) :: fr (nrmtmax), gr (nrmtmax), cf (3, nrmtmax)
      character(1024) :: message
! find the muffin-tin charges
      chgmttot = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               fr (ir) = rhomt (1, ir, ias) * spr (ir, is) ** 2
            End Do
            Call fderiv (-1, nrmt(is), spr(:, is), fr, gr, cf)
            chgmt (ias) = fourpi * y00 * gr (nrmt(is))
            chgmttot = chgmttot + chgmt (ias)
         End Do
      End Do
! find the interstitial charge
      sum = 0.d0
      Do ir = 1, ngrtot
         sum = sum + rhoir (ir) * cfunir (ir)
      End Do
      chgir = sum * omega / dble (ngrtot)
! total calculated charge
      chgcalc = chgmttot + chgir
      t1 = chgtot / chgcalc
      If (Abs(t1-1.d0) .Gt. input%groundstate%epschg) Then
         call warning('Warning(charge):')
         write(message,'(" Total charge density incorrect for s.c. loop ", I5)') iscl
         call warning(message)
         write(message,'(" Calculated : ", G18.10)') chgcalc
         call warning(message)
         write(message,'(" Required   : ", G18.10)') chgtot
         call warning(message)
      End If
      Return
End Subroutine
!EOC
