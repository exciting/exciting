!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmdedn (dedn)
! calculates derivative of total energy w.r.t. occupation numbers
      Use modinput
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (Out) :: dedn (nstsv, nkpt)
! allocatable
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: c (:, :)
      Integer :: ik, ist
      Allocate (evecsv(nstsv, nstsv), c(nstsv, nstsv))
      dedn (:, :) = 0.d0
      Do ik = 1, nkpt
! get evecsv from a file
         Call getevecsv (vkl(:, ik), evecsv)
! kinetic contribution
         Call zgemm ('C', 'N', nstsv, nstsv, nstsv, zone, evecsv, &
        & nstsv, dkdc(:, :, ik), nstsv, zzero, c, nstsv)
         Do ist = 1, nstsv
! include Coulomb contribution
            dedn (ist, ik) = dedn (ist, ik) - (dble(c(ist, &
           & ist))+dble(vclmat(ist, ist, ik)))
         End Do
      End Do
! add exchange correlation contribution
      Call rdmdexcdn (dedn)
! add entropic contribution if needed
      If (input%groundstate%RDMFT%rdmtemp .Gt. 0.d0) Then
         Call rdmdsdn (dedn)
      End If
      Deallocate (evecsv, c)
      Return
End Subroutine
