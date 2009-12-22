!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine rdmvaryn
! calculate new occupation numbers using derivatives of total energy
      Use modinput
      Use modmain
      Implicit None
! local variables
      Integer, Parameter :: maxit = 10000
      Integer :: it, ik, ist
      Real (8), Parameter :: eps = 1.d-12
      Real (8) :: tau, sum, gs, gsp, dgs
      Real (8) :: kapa, dkapa, t1
! allocatable arrays
      Real (8), Allocatable :: dedn (:, :)
      Real (8), Allocatable :: gamma (:, :)
      Allocate (dedn(nstsv, nkpt))
      Allocate (gamma(nstsv, nkpt))
! add constant to occupancies for charge conservation
      sum = 0.d0
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            sum = sum + wkpt (ik) * occsv (ist, ik)
         End Do
      End Do
      t1 = (chgval-sum) / dble (nstsv)
      occsv (:, :) = occsv (:, :) + t1
! redistribute charge so that occupancies are in the interval [0,occmax]
      sum = 0.d0
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            If (occsv(ist, ik) .Gt. occmax) Then
               sum = sum + wkpt (ik) * (occsv(ist, ik)-occmax)
               occsv (ist, ik) = occmax
            End If
            If (occsv(ist, ik) .Lt. 0.d0) Then
               sum = sum + wkpt (ik) * occsv (ist, ik)
               occsv (ist, ik) = 0.d0
            End If
         End Do
      End Do
      Do ist = 1, nstsv
         Do ik = 1, nkpt
            If (sum .Gt. 0.d0) Then
               t1 = wkpt (ik) * (occmax-occsv(ist, ik))
               t1 = Min (t1, sum)
               occsv (ist, ik) = occsv (ist, ik) + t1 / wkpt (ik)
               sum = sum - t1
            Else
               t1 = wkpt (ik) * occsv (ist, ik)
               t1 = Min (t1,-sum)
               occsv (ist, ik) = occsv (ist, ik) - t1 / wkpt (ik)
               sum = sum + t1
            End If
         End Do
      End Do
! get the derivatives
      Call rdmdedn (dedn)
! find suitable value of kapa such that sum of gamma is 0
      gsp = 0.d0
      kapa = 0.d0
      dkapa = 0.1d0
      Do it = 1, maxit
         gs = 0.d0
         sum = 0.d0
         Do ik = 1, nkpt
            Do ist = 1, nstsv
               t1 = dedn (ist, ik) - kapa
               If (t1 .Gt. 0.d0) Then
                  gamma (ist, ik) = t1 * (occmax-occsv(ist, ik))
               Else
                  gamma (ist, ik) = t1 * occsv (ist, ik)
               End If
               gs = gs + wkpt (ik) * gamma (ist, ik)
               sum = sum + wkpt (ik) * gamma (ist, ik) ** 2
            End Do
         End Do
         sum = Sqrt (sum)
         sum = Max (sum, 1.d0)
         t1 = Abs (gs) / sum
         If (t1 .Lt. eps) Go To 10
         If (it .Ge. 2) Then
            dgs = gs - gsp
            If (gs*dgs .Gt. 0.d0) dkapa = - dkapa
            If (gs*gsp .Lt. 0.d0) Then
               dkapa = 0.5d0 * dkapa
            Else
               dkapa = 1.1d0 * dkapa
            End If
         End If
         gsp = gs
         kapa = kapa + dkapa
      End Do
      Write (*,*)
      Write (*, '("Error(rdmdedn): could not find offset")')
      Write (*,*)
      Stop
10    Continue
! normalize gamma if sum of squares is greater than 1
      sum = 0.d0
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            sum = sum + wkpt (ik) * gamma (ist, ik) ** 2
         End Do
      End Do
      If (sum .Gt. 1.d0) Then
         t1 = 1.d0 / Sqrt (sum)
         gamma (:, :) = t1 * gamma (:, :)
      End If
! find step size which keeps occupancies in the interval [0,occmax]
      tau = input%groundstate%RDMFT%taurdmn
20    Continue
      If (Abs(tau) .Lt. eps) Go To 30
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            t1 = occsv (ist, ik) + tau * gamma (ist, ik)
            If (gamma(ist, ik) .Gt. 0.d0) Then
               If (t1 .Gt. occmax+eps) Then
                  tau = 0.75d0 * tau
                  Go To 20
               End If
            End If
            If (gamma(ist, ik) .Lt. 0.d0) Then
               If (t1 .Lt.-eps) Then
                  tau = 0.75d0 * tau
                  Go To 20
               End If
            End If
         End Do
      End Do
30    Continue
! update occupancies
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            occsv (ist, ik) = occsv (ist, ik) + tau * gamma (ist, ik)
         End Do
      End Do
! write derivatives and occupancies to a file
      Call rdmwritededn (dedn)
      Deallocate (dedn, gamma)
      Return
End Subroutine
