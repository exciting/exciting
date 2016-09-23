!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rfmtctof
! !INTERFACE:
!
!
Subroutine rfmtctof (rfmt)
! !INPUT/OUTPUT PARAMETERS:
      Use modinput
!   rfmt : real muffin-tin function (in,real(lmmaxvr,nrmtmax,natmtot))
! !DESCRIPTION:
!   Converts a real muffin-tin function from a coarse to a fine radial mesh by
!   using cubic spline interpolation. See routines {\tt rfinterp} and
!   {\tt spline}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (Inout) :: rfmt (lmmaxvr, nrmtmax, natmtot)
! local variables
      Integer :: is, ia, ias, ld, lm
      ld = lmmaxvr * input%groundstate%lradstep
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! interpolate with a clamped spline
            Do lm = 1, lmmaxvr
               Call rfinterp (nrcmt(is), rcmt(:, is), ld, rfmt(lm, 1, &
              & ias), nrmt(is), spr(:, is), lmmaxvr, rfmt(lm, 1, ias))
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
