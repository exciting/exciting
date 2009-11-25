!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmdkdc
! calculates the derivative of kinetic energy w.r.t. evecsv
      Use modmain
      Implicit None
! allocatable arrays
      Complex (8), Allocatable :: evecsv (:, :)
      Integer :: ik
      Allocate (evecsv(nstsv, nstsv))
      Do ik = 1, nkpt
         Call getevecsv (vkl(:, ik), evecsv)
         Call zgemm ('N', 'N', nstsv, nstsv, nstsv, zone, kinmatc(:, :, &
        & ik), nstsv, evecsv, nstsv, zzero, dkdc(:, :, ik), nstsv)
      End Do
      Deallocate (evecsv)
      Return
End Subroutine
