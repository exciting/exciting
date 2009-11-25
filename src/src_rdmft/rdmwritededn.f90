!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine rdmwritededn (dedn)
! writes derivative of total energy w.r.t. occupancies to file
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: dedn (nstsv, nkpt)
! local variables
      Integer :: ik, ist
      Open (50, File='RDM_DEDN.OUT', Action='WRITE', Form='FORMATTED')
      Write (50, '(I6, " : nkpt")') nkpt
      Write (50, '(I6, " : nstsv")') nstsv
      Do ik = 1, nkpt
         Write (50,*)
         Write (50, '(I6, 3G18.10, " : k-point, vkl")') ik, vkl (:, ik)
         Write (50, '("     (state, occupancy and derivative below)")')
         Do ist = 1, nstsv
            Write (50, '(I6, 4G18.10)') ist, occsv (ist, ik), - dedn &
           & (ist, ik)
         End Do
      End Do
      Close (50)
      Return
End Subroutine
