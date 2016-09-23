!
!
!
! Copyright (C) 2008 T. Baldsiefen, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine rdmentropy
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: ik, ist
      Real (8) :: t1
      rdmentrpy = 0.d0
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            t1 = Max (occsv(ist, ik), input%groundstate%epsocc)
            t1 = Min (t1, occmax-input%groundstate%epsocc)
            rdmentrpy = rdmentrpy - wkpt (ik) * &
           & (t1*Log(t1/occmax)+(occmax-t1)*Log(1.d0-t1/occmax))
         End Do
      End Do
      rdmentrpy = kboltz * rdmentrpy
      Return
End Subroutine
