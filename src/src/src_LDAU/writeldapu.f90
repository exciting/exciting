!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine writeldapu
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ia, ias, ispn, jspn
      Integer :: l, m1, m2, lm1, lm2
      Open (50, File='DMATLU'//trim(filext), Action='WRITE', Form='FORM&
     &ATTED')
      Do is = 1, nspecies
         l = llu (is)
         If (l .Ge. 0) Then
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Write (50,*)
               Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') is, trim &
              & (input%structure%speciesarray(is)%species%chemicalSymbol), ia
               Write (50, '(" l = ", I2)') l
               Do ispn = 1, nspinor
                  Do jspn = 1, nspinor
                     Write (50, '(" ispn = ", I1, ", jspn = ", I1)') &
                    & ispn, jspn
                     Do m1 = - l, l
                        lm1 = idxlm (l, m1)
                        Do m2 = - l, l
                           lm2 = idxlm (l, m2)
                           Write (50, '("  m1 = ", I2, ", m2 = ", I2, "&
                          & : ", 2G18.10)') m1, m2, dmatlu (lm1, lm2, &
                          & ispn, jspn, ias)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End If
      End Do
      Close (50)
      Open (50, File='VMATLU'//trim(filext), Action='WRITE', Form='FORM&
     &ATTED')
      Do is = 1, nspecies
         l = llu (is)
         If (l .Ge. 0) Then
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Write (50,*)
               Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') is, trim &
              & (input%structure%speciesarray(is)%species%chemicalSymbol), ia
               Write (50, '(" l = ", I2)') l
               Do ispn = 1, nspinor
                  Do jspn = 1, nspinor
                     Write (50, '(" ispn = ", I1, ", jspn = ", I1)') &
                    & ispn, jspn
                     Do m1 = - l, l
                        lm1 = idxlm (l, m1)
                        Do m2 = - l, l
                           lm2 = idxlm (l, m2)
                           Write (50, '("  m1 = ", I2, ", m2 = ", I2, "&
                          & : ", 2G18.10)') m1, m2, vmatlu (lm1, lm2, &
                          & ispn, jspn, ias)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End If
      End Do
      Close (50)
      If (ldapu .Eq. 3) Then
         Open (50, File='ALPHALU'//trim(filext), Action='WRITE', Form='&
        &FORMATTED')
         Do is = 1, nspecies
            l = llu (is)
            If (l .Ge. 0) Then
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  Write (50,*)
                  Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') is, trim &
                 & (input%structure%speciesarray(is)%species%chemicalSymbol), ia
                  Write (50, '(" l = ", I2)') l
                  Write (50, '(" alpha = ", G18.10)') alphalu (ias)
               End Do
            End If
         End Do
         Close (50)
      End If
      Return
End Subroutine
