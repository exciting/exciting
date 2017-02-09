
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: copyfilesq0
! !INTERFACE:
Subroutine copyfilesq0
      Use modinput
! !USES:
      Use modmain
      Use modxs
      Use m_getapwcmt
      Use m_getlocmt
! !DESCRIPTION:
!   If a finite momentum transfer Q-point is the Gamma-point, eigenvalues
!   eigenvectors, occupancies and muffin-tin expansion coefficients are
!   identical to those corresponding to the unshifted mesh. Files are copied
!   using the associated generic routines for ISO compatibility, or links
!   are generated if ISO compatibility is dropped.
!
! !REVISION HISTORY:
!   Created January 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Integer, Parameter :: iq = 1
      Integer :: ik
      Complex (8), Allocatable :: evecfvt (:, :, :)
      Complex (8), Allocatable :: apwlm (:, :, :, :), lolm (:, :, :, :)
      Call init1offs (qvkloff(1, iq))
      Allocate (evecfvt(nmatmax, nstfv, nspnfv), evecsv(nstsv, nstsv))
      Allocate (apwlm(nstfv, apwordmax, lmmaxapw, natmtot))
      Allocate (lolm(nstfv, nlomax,-lolmax:lolmax, natmtot))
      Do ik = 1, nkpt
     ! read files
        If (skipgnd) Then
           filext = '.OUT'
        Else
           filext = '_QMT001.OUT'  
        End If  
         Call getevecfv (vkl(1, ik), vgkl(1, 1, 1, ik), evecfvt)
         Call getevecsv (vkl(1, ik), evecsv)
         Call getevalsv (vkl(1, ik), evalsv(1, ik))
         Call getoccsv (vkl(1, ik), occsv(1, ik))
         Call getapwcmt (iq, ik, 1, nstfv, input%groundstate%lmaxapw, &
        & apwlm)
         Call getlocmt (iq, ik, 1, nstfv, lolm)
     ! write files
         filext = '_QMT000.OUT'
         Call putevecfv (ik, evecfvt)
         Call putevecsv (ik, evecsv)
         Call putevalsv (ik, evalsv(1, ik))
         Call putoccsv (ik, occsv(1, ik))
         Call putapwcmt ('APWCMT_QMT000.OUT', ik, vkl(:, ik), vql(:, &
        & iq), apwlm)
         Call putlocmt ('LOCMT_QMT000.OUT', ik, vkl(:, ik), vql(:, iq), &
        & lolm)
      End Do
  ! read files
      If (skipgnd) Then
           filext = '.OUT'
      Else
           filext = '_QMT001.OUT'  
      End If  
      Call readfermi
  ! write files
      filext = '_QMT000.OUT'
      Call writeeval
      Call writefermi
      Deallocate (evecfvt)
      Deallocate (evecsv)
      Deallocate (apwlm)
      Deallocate (lolm)
End Subroutine copyfilesq0
!EOC
