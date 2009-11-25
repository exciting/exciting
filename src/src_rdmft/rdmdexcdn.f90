!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmdexcdn (dedn)
! calculates derivative of exchange-correlation energy w.r.t. occupation numbers
      Use modinput
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (Inout) :: dedn (nstsv, nkpt)
! local variables
      Integer :: ik1, ik2, ik3
      Integer :: ist1, ist2, iv (3)
! parameter for calculating the functional derivatives
      Real (8), Parameter :: eps = 1.d-12
      Real (8) :: t1, t2, t3, t4
! external functions
      Real (8) :: r3taxi
      External r3taxi
      If (input%groundstate%RDMFT%rdmxctype .Eq. 0) Return
! calculate the pre-factor
      If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
         t1 = 1.d0 / occmax
      Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) Then
         If (associated(input%groundstate%spin)) Then
            t1 = input%groundstate%RDMFT%rdmalpha
         Else
            t1 = 2.d0 * input%groundstate%RDMFT%rdmalpha * (0.25d0) ** &
           & input%groundstate%RDMFT%rdmalpha
         End If
      Else
         Write (*,*)
         Write (*, '("Error(rdmdexcdn): rdmxctype not defined : ", I8)') input%groundstate%RDMFT%rdmxctype
         Write (*,*)
      End If
      Do ik1 = 1, nkpt
         Do ist1 = 1, nstsv
            Do ik2 = 1, nkptnr
! find the equivalent reduced k-point
               iv (:) = ivknr (:, ik2)
               ik3 = ikmap (iv(1), iv(2), iv(3))
               Do ist2 = 1, nstsv
! Hartree-Fock functional
                  If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
                     t2 = t1 * occsv (ist2, ik3)
! SDLG functional
                  Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) &
                 & Then
                     If ((ist1 .Eq. ist2) .And. (r3taxi(vkl(1, ik1), &
                    & vklnr(1, ik2)) .Lt. input%structure%epslat)) Then
                        t2 = (1.d0/occmax) * occsv (ist2, ik3)
                     Else
                        t3 = Max (occsv(ist1, ik1), eps)
                        t4 = Max (occsv(ist2, ik3), eps)
                        t2 = t1 * &
                       & (t4**input%groundstate%RDMFT%rdmalpha) / &
                       & (t3**(1.d0-input%groundstate%RDMFT%rdmalpha))
                     End If
                  End If
                  dedn (ist1, ik1) = dedn (ist1, ik1) + t2 * vnlrdm &
                 & (ist1, ik1, ist2, ik2)
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
