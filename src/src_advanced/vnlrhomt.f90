!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: vnlrhomt
! !INTERFACE:
!
!
Subroutine vnlrhomt (tsh, is, wfmt1, wfmt2, zrhomt)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh    : .true. if the density is to be in spherical harmonics (in,logical)
!   is     : species number (in,integer)
!   wfmt1  : muffin-tin part of wavefunction 1 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfmt2  : muffin-tin part of wavefunction 2 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   zrhomt : muffin-tin charge density in spherical harmonics/coordinates
!            (out,complex(lmmaxvr,nrcmtmax))
! !DESCRIPTION:
!   Calculates the complex overlap density in a single muffin-tin from two input
!   wavefunctions expressed in spherical coordinates. If {\tt tsh} is
!   {\tt .true.} then the output density is converted to a spherical harmonic
!   expansion. See routine {\tt vnlrho}.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tsh
      Integer, Intent (In) :: is
      Complex (8), Intent (In) :: wfmt1 (lmmaxvr, nrcmtmax)
      Complex (8), Intent (In) :: wfmt2 (lmmaxvr, nrcmtmax)
      Complex (8), Intent (Out) :: zrhomt (lmmaxvr, nrcmtmax)
! local variables
      Integer :: irc
! allocatable arrays
      Complex (8), Allocatable :: zfmt (:, :)
      If (tsh) Then
! output density in spherical harmonics
         Allocate (zfmt(lmmaxvr, nrcmtmax))
         Do irc = 1, nrcmt (is)
            zfmt (:, irc) = conjg (wfmt1(:, irc)) * wfmt2 (:, irc)
         End Do
         Call zgemm ('N', 'N', lmmaxvr, nrcmt(is), lmmaxvr, zone, &
        & zfshtvr, lmmaxvr, zfmt, lmmaxvr, zzero, zrhomt, lmmaxvr)
         Deallocate (zfmt)
      Else
! output density in spherical coordinates
         Do irc = 1, nrcmt (is)
            zrhomt (:, irc) = conjg (wfmt1(:, irc)) * wfmt2 (:, irc)
         End Do
      End If
      Return
End Subroutine
!EOC
