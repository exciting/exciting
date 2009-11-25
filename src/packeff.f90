!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: packeff
! !INTERFACE:
!
!
Subroutine packeff (tpack, n, nu)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tpack : .true. for packing, .false. for unpacking (in,logical)
!   n     : total number of real values stored (out,integer)
!   nu    : packed potential (inout,real(*))
! !DESCRIPTION:
!   Packs/unpacks the muffin-tin and interstitial parts of the effective
!   potential and magnetic field into/from the single array {\tt nu}. This array
!   can then be passed directly to the mixing routine. See routine {\tt rfpack}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tpack
      Integer, Intent (Out) :: n
      Real (8), Intent (Inout) :: nu (*)
! local variables
      Integer :: idm, ias, lm1, lm2
      Integer :: ispn, jspn
      n = 0
      Call rfpack (tpack, n, 1, veffmt, veffir, nu)
      Do idm = 1, ndmag
         Call rfpack (tpack, n, 1, bxcmt(:, :, :, idm), bxcir(:, idm), &
        & nu)
      End Do
! pack the LDA+U potential if required
      If (ldapu .Ne. 0) Then
         Do ias = 1, natmtot
            Do ispn = 1, nspinor
               Do jspn = 1, nspinor
                  Do lm1 = 1, lmmaxlu
                     Do lm2 = 1, lmmaxlu
                        n = n + 1
                        If (tpack) Then
                           nu (n) = dble (vmatlu(lm1, lm2, ispn, jspn, &
                          & ias))
                           n = n + 1
                           nu (n) = aimag (vmatlu(lm1, lm2, ispn, jspn, &
                          & ias))
                        Else
                           vmatlu (lm1, lm2, ispn, jspn, ias) = cmplx &
                          & (nu(n), nu(n+1), 8)
                           n = n + 1
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If
      Return
End Subroutine
!EOC
