!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genrmesh
! !INTERFACE:
!
!
Subroutine genrmesh
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Generates the coarse and fine radial meshes for each atomic species in the
!   crystal. Also determines which points are in the inner part of the
!   muffin-tin using the value of {\tt radfinr}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ir, irc
      Real (8) :: t1, t2
! estimate the number of radial mesh points to infinity
      spnrmax = 1
      Do is = 1, nspecies
! logarithmic mesh
         t1 = Log (sprmax(is)/sprmin(is)) / Log (rmt(is)/sprmin(is))
         t2 = dble (nrmt(is)-1) * t1
         spnr (is) = Nint (t2) + 1
         spnrmax = Max (spnrmax, spnr(is))
      End Do
! allocate the global radial mesh arrays
      If (allocated(spr)) deallocate (spr)
      Allocate (spr(spnrmax, nspecies))
      If (allocated(rcmt)) deallocate (rcmt)
      Allocate (rcmt(nrcmtmax, nspecies))
! generate the radial meshes
      Do is = 1, nspecies
         t1 = 1.d0 / dble (nrmt(is)-1)
! logarithmic mesh
         t2 = Log (rmt(is)/sprmin(is))
         Do ir = 1, spnr (is)
            spr (ir, is) = sprmin (is) * Exp (dble(ir-1)*t1*t2)
         End Do
      End Do
! find the inner part of the muffin-tin (where rho is calculated with lmaxinr)
      Do is = 1, nspecies
         t1 = input%groundstate%fracinr * rmt (is)
         nrmtinr (is) = nrmt (is)
         Do ir = 1, nrmt (is)
            If (spr(ir, is) .Gt. t1) Then
               nrmtinr (is) = ir
               Go To 10
            End If
         End Do
10       Continue
      End Do
! set up the coarse radial meshes
      Do is = 1, nspecies
         irc = 0
         Do ir = 1, nrmt (is), input%groundstate%lradstep
            irc = irc + 1
            rcmt (irc, is) = spr (ir, is)
         End Do
      End Do
      Return
End Subroutine
!EOC
