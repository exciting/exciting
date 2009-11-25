!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writeeval
! !INTERFACE:
!
!
Subroutine writeeval
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Outputs the second-variational eigenvalues and occupation numbers to the
!   file {\tt EIGVAL.OUT}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, ist, is, ia, ias
! write out the valence eigenvalues
      Open (50, File='EIGVAL'//trim(filext), Action='WRITE', Form='FORM&
     &ATTED')
      Write (50, '(I6, " : nkpt")') nkpt
      Write (50, '(I6, " : nstsv")') nstsv
      Do ik = 1, nkpt
         Write (50,*)
         Write (50, '(I6, 3G18.10, " : k-point, vkl")') ik, vkl (:, ik)
         Write (50, '(" (state, eigenvalue and occupancy below)")')
         Do ist = 1, nstsv
            Write (50, '(I6, 2G18.10)') ist, evalsv (ist, ik), occsv &
           & (ist, ik)
         End Do
         Write (50,*)
      End Do
      Close (50)
! write out the core eigenvalues
      Open (50, File='EVALCORE'//trim(filext), Action='WRITE', Form='FO&
     &RMATTED')
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Write (50,*)
            Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') &
           & is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol), &
           & ia
            Do ist = 1, spnst (is)
               If (spcore(ist, is)) Then
                  Write (50, '(" n = ", I2, ", l = ", I2, ", k = ", I2,&
                 & " : ", G18.10)') spn (ist, is), spl (ist, is), spk &
                 & (ist, is), evalcr (ist, ias)
               End If
            End Do
         End Do
      End Do
      Close (50)
      Return
End Subroutine
!EOC
