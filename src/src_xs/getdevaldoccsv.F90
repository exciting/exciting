!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getdevaldoccsv (iq, ik, ikq, l1, u1, l2, u2, devalsv, &
& doccsv, scissv)
  ! xssave0 has to be called in advance.
      Use modinput
      Use modmain
      Use modxs
      Use m_genfilname
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik, ikq, l1, u1, l2, u2
      Real (8), Intent (Out) :: devalsv (u1-l1+1, u2-l2+1)
      Real (8), Intent (Out) :: doccsv (u1-l1+1, u2-l2+1)
      Real (8), Intent (Out) :: scissv (u1-l1+1, u2-l2+1)
  ! local variables
      Integer :: ist, jst, iqt
      Real (8), Allocatable :: e0 (:), e (:), o0 (:), o (:)
      iqt = iq
      Allocate (e0(nstsv), e(nstsv), o0(nstsv), o(nstsv))
  ! eigenvalues and occupancies for k+q-point
      Call getevalsv (vkl(1, ikq), e)
      Call getoccsv (vkl(1, ikq), o)
  ! eigenvalues and occupancies for k-point
      Call getevalsv0 (vkl0(1, ik), e0)
      Call getoccsv0 (vkl0(1, ik), o0)
  ! scissors correction
      scissv (:, :) = 0.d0
      Do ist = l1, u1
         Do jst = l2, u2
        ! neglect contributions above cutoff
            If ((e0(ist) .Lt. input%xs%emaxdf) .And. (e(jst) .Lt. &
           & input%xs%emaxdf)) Then
               devalsv (ist-l1+1, jst-l2+1) = e0 (ist) - e (jst)
               doccsv (ist-l1+1, jst-l2+1) = o0 (ist) - o (jst)
               If ((e0(ist) .Le. efermi) .And. (e(jst) .Gt. efermi)) &
              & scissv (ist-l1+1, jst-l2+1) = - input%xs%scissor
               If ((e0(ist) .Gt. efermi) .And. (e(jst) .Le. efermi)) &
              & scissv (ist-l1+1, jst-l2+1) = input%xs%scissor
            Else
               If (input%xs%dbglev .Gt. 1) Then
                  Write (*, '("Info(getdevaldoccsv): cutoff applied: iq&
                 &, ik, ist, jst, energies:", 4i6, 2g18.10)') iq, ik, &
                 & ist, jst, e0 (ist), e (jst)
               End If
           ! set energy difference to arbitrary number
               devalsv (ist-l1+1, jst-l2+1) = 1.d0
           ! set occupation number difference to zero as a factor in the
           ! terms of the dielectric function
               doccsv (ist-l1+1, jst-l2+1) = 0.d0
            End If
         End Do
      End Do
      Deallocate (e0, e, o0, o)
End Subroutine getdevaldoccsv
