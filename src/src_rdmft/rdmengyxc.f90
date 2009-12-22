!
!
!
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine rdmengyxc
! calculate the RDMFT exchange-correlation energy
      Use modinput
      Use modmain
      Implicit None
! local variables
      Integer :: ik1, ik2, ik3
      Integer :: ist1, ist2, iv (3)
      Real (8) :: t1, t2, t3
! external functions
      Real (8) :: r3taxi, rfmtinp
      External r3taxi, rfmtinp
! calculate the prefactor
      If (input%groundstate%RDMFT%rdmxctype .Eq. 0) Then
         engyx = 0.d0
         Return
      Else If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
         t1 = 0.5d0 / occmax
      Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) Then
         If (associated(input%groundstate%spin)) Then
            t1 = 0.5d0
         Else
            t1 = (0.25d0) ** input%groundstate%RDMFT%rdmalpha
         End If
      Else
         Write (*,*)
         Write (*, '("Error(rdmengyxc): rdmxctype not defined : ", I8)') input%groundstate%RDMFT%rdmxctype
         Write (*,*)
         Stop
      End If
! exchange-correlation energy
      engyx = 0.d0
      Do ik1 = 1, nkpt
         Do ist1 = 1, nstsv
            Do ik2 = 1, nkptnr
! find the equivalent reduced k-point
               iv (:) = ivknr (:, ik2)
               ik3 = ikmap (iv(1), iv(2), iv(3))
               Do ist2 = 1, nstsv
! Hartree-Fock functional
                  If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
                     t2 = t1 * wkpt (ik1) * occsv (ist1, ik1) * occsv &
                    & (ist2, ik3)
! SDLG functional
                  Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) &
                 & Then
                     t3 = occsv (ist1, ik1) * occsv (ist2, ik3)
                     If ((ist1 .Eq. ist2) .And. (r3taxi(vkl(1, ik1), &
                    & vklnr(1, ik2)) .Lt. input%structure%epslat)) Then
                        t2 = (0.5d0/occmax) * wkpt (ik1) * t3
                     Else
                        If (t3 .Gt. 0.d0) Then
                           t2 = t1 * wkpt (ik1) * t3 ** &
                          & input%groundstate%RDMFT%rdmalpha
                        Else
                           t2 = 0.d0
                        End If
                     End If
                  End If
                  engyx = engyx - t2 * vnlrdm (ist1, ik1, ist2, ik2)
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
