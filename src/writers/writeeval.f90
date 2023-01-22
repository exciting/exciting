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
      Use Fox_wxml
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
      Type (xmlf_t), Save :: xf
! write out the valence eigenvalues
      Open (50, File='EIGVAL'//trim(filext), Action='WRITE', Form='FORM&
     &ATTED')
      Write (50, '(I6, " : nkpt")') nkpt
      Write (50, '(I6, " : nstsv")') nstsv
      Call xml_OpenFile("eigval.xml", xf, replace = .True., pretty_print = .True.)
      Call xml_NewElement(xf, "eigval")
      Do ik = 1, nkpt
         Write (50,*)
         Write (50, '(I6, 3G18.10, " : k-point, vkl")') ik, vkl (:, ik)
         Write (50, '(" (state, eigenvalue and occupancy below)")')
         Call xml_NewElement(xf, "kpt")
         Call xml_AddAttribute(xf, "ik", ik)
         Call xml_AddAttribute(xf, "vkl", vkl(:, ik))
         Do ist = 1, nstsv
            Write (50, '(I6, 2G18.10)') ist, evalsv (ist, ik), occsv &
           & (ist, ik)
           Call xml_NewElement(xf, "state")
           Call xml_AddAttribute(xf, "ist", ist)
           Call xml_AddAttribute(xf, "eigenvalue", evalsv (ist, ik))
           Call xml_AddAttribute(xf, "occupancy", occsv(ist, ik))
           Call xml_EndElement(xf, "state")
         End Do
         Write (50,*)
         Call xml_EndElement(xf, "kpt")
      End Do
      Close (50)
      Call xml_EndElement(xf, "eigval")
      Call xml_Close(xf)
! write out the core eigenvalues
      Open (50, File='EVALCORE'//trim(filext), Action='WRITE', Form='FO&
     &RMATTED')
      Call xml_OpenFile("evalcore.xml", xf, replace = .True., pretty_print = .True.)
      Call xml_NewElement(xf, "evalcore")
      Do is = 1, nspecies
         Call xml_NewElement(xf, "species")
         Call xml_AddAttribute(xf, "is", is)
         Call xml_AddAttribute(xf, "chemicalSymbol", &
            &trim(input%structure%speciesarray(is)%species%chemicalSymbol))
         Do ia = 1, natoms (is)
            Call xml_NewElement(xf, "atom")
            Call xml_AddAttribute(xf, "ia", ia)
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
                 Call xml_NewElement(xf, "state")
                 Call xml_AddAttribute(xf, "ist", ist)
                 Call xml_AddAttribute(xf, "n", spn(ist, is))
                 Call xml_AddAttribute(xf, "l", spl(ist, is))
                 Call xml_AddAttribute(xf, "k", spl(ist, is))
                 Call xml_AddAttribute(xf, "eigenvalue", evalcr (ist, ias))
                 Call xml_EndElement(xf, "state")
               End If
            End Do !ist
            Call xml_EndElement(xf, "atom")
         End Do !ia
         Call xml_EndElement(xf, "species")
      End Do !is
      Close (50)
      Call xml_EndElement(xf, "evalcore")
      Call xml_Close(xf)
      Return
End Subroutine
!EOC
