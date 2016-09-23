!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine oepmagmt (tsh, is, wfmt1, wfmt2, wfmt3, wfmt4, zvfmt)
      Use modmain
      Implicit None
! arguments
      Logical, Intent (In) :: tsh
      Integer, Intent (In) :: is
      Complex (8), Intent (In) :: wfmt1 (lmmaxvr, nrcmtmax)
      Complex (8), Intent (In) :: wfmt2 (lmmaxvr, nrcmtmax)
      Complex (8), Intent (In) :: wfmt3 (lmmaxvr, nrcmtmax)
      Complex (8), Intent (In) :: wfmt4 (lmmaxvr, nrcmtmax)
      Complex (8), Intent (Out) :: zvfmt (lmmaxvr, nrcmtmax, ndmag)
! local variables
      Integer :: nrc
! allocatable arrays
      Complex (8), Allocatable :: zfmt (:, :, :)
      Allocate (zfmt(lmmaxvr, nrcmtmax, 2))
! muffin-tin part
      nrc = nrcmt (is)
! up-up spin density
      Call vnlrhomt (tsh, is, wfmt1, wfmt3, zfmt(:, :, 1))
! dn-dn spin density
      Call vnlrhomt (tsh, is, wfmt2, wfmt4, zfmt(:, :, 2))
! calculate the z-component of mangetisation: up-up - dn-dn
      zvfmt (:, 1:nrc, ndmag) = zfmt (:, 1:nrc, 1) - zfmt (:, 1:nrc, 2)
! non-collinear case
      If (ncmag) Then
! up-dn spin density
         Call vnlrhomt (tsh, is, wfmt1, wfmt4, zfmt(:, :, 1))
! dn-up spin density
         Call vnlrhomt (tsh, is, wfmt2, wfmt3, zfmt(:, :, 2))
! calculate the x-component: up-dn + dn-up
         zvfmt (:, 1:nrc, 1) = zfmt (:, 1:nrc, 1) + zfmt (:, 1:nrc, 2)
! calculate the y-component: i*(dn-up - up-dn)
         zvfmt (:, 1:nrc, 2) = zi * (zfmt(:, 1:nrc, 2)-zfmt(:, 1:nrc, &
        & 1))
      End If
      Deallocate (zfmt)
      Return
End Subroutine
