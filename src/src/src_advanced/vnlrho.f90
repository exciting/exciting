!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: vnlrho
! !INTERFACE:
!
!
Subroutine vnlrho (tsh, wfmt1, wfmt2, wfir1, wfir2, zrhomt, zrhoir)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh    : .true. if the muffin-tin density is to be in spherical harmonics
!            (in,logical)
!   wfmt1  : muffin-tin part of wavefunction 1 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfmt2  : muffin-tin part of wavefunction 2 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfir1  : interstitial wavefunction 1 (in,complex(ngrtot))
!   wfir2  : interstitial wavefunction 2 (in,complex(ngrtot))
!   zrhomt : muffin-tin charge density in spherical harmonics/coordinates
!            (out,complex(lmmaxvr,nrcmtmax,natmtot))
!   zrhoir : interstitial charge density (out,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the complex overlap charge density from two input wavefunctions:
!   $$ \rho({\bf r})\equiv\Psi_1^{\dag}({\bf r})\cdot\Psi_2({\bf r}). $$
!   Note that the muffin-tin wavefunctions are provided in spherical coordinates
!   and the returned density is either in terms of spherical harmonic
!   coefficients or spherical coordinates when {\tt tsh} is {\tt .true.} or
!   {\tt .false.}, respectively. See also the routine {\tt vnlrhomt}.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tsh
      Complex (8), Intent (In) :: wfmt1 (lmmaxvr, nrcmtmax, natmtot, &
     & nspinor)
      Complex (8), Intent (In) :: wfmt2 (lmmaxvr, nrcmtmax, natmtot, &
     & nspinor)
      Complex (8), Intent (In) :: wfir1 (ngrtot, nspinor)
      Complex (8), Intent (In) :: wfir2 (ngrtot, nspinor)
      Complex (8), Intent (Out) :: zrhomt (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (Out) :: zrhoir (ngrtot)
! local variables
      Integer :: is, ia, ias, nrc, ir
! allocatable arrays
      Complex (8), Allocatable :: zfmt (:, :)
      If (associated(input%groundstate%spin)) allocate (zfmt(lmmaxvr, &
     & nrcmtmax))
! muffin-tin part
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call vnlrhomt (tsh, is, wfmt1(:, :, ias, 1), wfmt2(:, :, &
           & ias, 1), zrhomt(:, :, ias))
            If (associated(input%groundstate%spin)) Then
! spin-polarised
               Call vnlrhomt (tsh, is, wfmt1(:, :, ias, 2), wfmt2(:, :, &
              & ias, 2), zfmt)
               zrhomt (:, 1:nrc, ias) = zrhomt (:, 1:nrc, ias) + zfmt &
              & (:, 1:nrc)
            End If
         End Do
      End Do
! interstitial part
      If (associated(input%groundstate%spin)) Then
! spin-polarised
         Do ir = 1, ngrtot
            zrhoir (ir) = conjg (wfir1(ir, 1)) * wfir2 (ir, 1) + conjg &
           & (wfir1(ir, 2)) * wfir2 (ir, 2)
         End Do
      Else
! spin-unpolarised
         Do ir = 1, ngrtot
            zrhoir (ir) = conjg (wfir1(ir, 1)) * wfir2 (ir, 1)
         End Do
      End If
      If (associated(input%groundstate%spin)) deallocate (zfmt)
      Return
End Subroutine
!EOC
