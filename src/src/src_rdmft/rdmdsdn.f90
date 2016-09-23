!
!
!
! Copyright (C) 2008 T. Baldsiefen, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine rdmdsdn (dedn)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (Inout) :: dedn (nstsv, nkpt)
! local variables
      Integer :: ik, ist
      Real (8) :: t1
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            t1 = Max (occsv(ist, ik), input%groundstate%epsocc)
            t1 = Min (t1, occmax-input%groundstate%epsocc)
            dedn (ist, ik) = dedn (ist, ik) - &
           & input%groundstate%RDMFT%rdmtemp * kboltz * Log &
           & (t1/(occmax-t1))
         End Do
      End Do
      Return
End Subroutine
