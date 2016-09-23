!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmdedc (ik, dedc)
! calculate the derivative of total energy w.r.t. evecsv
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (Out) :: dedc (nstsv, nstsv)
! local variables
      Integer :: ist
! allocatable arrays
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: c (:, :)
! allocate local arrays
      Allocate (evecsv(nstsv, nstsv))
      Allocate (c(nstsv, nstsv))
! get the eigenvectors from file
      Call getevecsv (vkl(:, ik), evecsv)
! kinetic and Coulomb potential contribution
      Call zgemm ('N', 'N', nstsv, nstsv, nstsv, zone, evecsv, nstsv, &
     & vclmat(:, :, ik), nstsv, zzero, c, nstsv)
      Do ist = 1, nstsv
         dedc (:, ist) = occsv (ist, ik) * (dkdc(:, ist, ik)+c(:, ist))
      End Do
! exchange-correlation contribution
      Call rdmdexcdc (ik, evecsv, dedc)
      Deallocate (evecsv, c)
      Return
End Subroutine
