!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmvaryc (sum)
! calculates new evecsv using the gradient of energy w.r.t. evecsv
      Use modinput
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (Out) :: sum
! local variables
      Integer :: ik, ist1, ist2
      Real (8) :: t1
      Complex (8) zt1
! allocatable arrays
      Complex (8), Allocatable :: dedc (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: evecsvt (:)
! external functions
      Real (8) :: dznrm2
      Complex (8) zdotc
      External dznrm2, zdotc
      Allocate (dedc(nstsv, nstsv, nkpt))
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
      Do ik = 1, nkpt
!$OMP CRITICAL
         Write (*, '("Info(rdmvaryc): ", I6, " of ", I6, " k-points")') &
        & ik, nkpt
!$OMP END CRITICAL
! compute the derivative w.r.t. evecsv
         Call rdmdedc (ik, dedc(:, :, ik))
      End Do
!$OMP END DO
!$OMP END PARALLEL
      Allocate (evecsv(nstsv, nstsv))
      Allocate (evecsvt(nstsv))
      sum = 0.d0
      Do ik = 1, nkpt
! get the eigenvectors from file
         Call getevecsv (vkl(:, ik), evecsv)
! calculate new evecsv
         evecsv (:, :) = evecsv (:, :) - &
        & input%groundstate%RDMFT%taurdmc * dedc (:, :, ik)
! othogonalise evecsv (Gram-Schmidt)
         Do ist1 = 1, nstsv
            evecsvt (:) = evecsv (:, ist1)
            Do ist2 = 1, ist1 - 1
               zt1 = zdotc (nstsv, evecsv(:, ist2), 1, evecsv(:, ist1), &
              & 1)
               evecsvt (:) = evecsvt (:) - zt1 * evecsv (:, ist2)
            End Do
            t1 = dznrm2 (nstsv, evecsvt, 1)
            t1 = 1.d0 / t1
            evecsv (:, ist1) = t1 * evecsvt (:)
         End Do
! write new evecsv to file
         Call putevecsv (ik, evecsv)
! convergence check
         Do ist1 = 1, nstsv
            Do ist2 = 1, nstsv
               zt1 = dedc (ist1, ist2, ik)
               sum = sum + wkpt (ik) * (dble(zt1)**2+aimag(zt1)**2)
            End Do
         End Do
! end loop over k-points
      End Do
      Deallocate (dedc, evecsv, evecsvt)
      Return
End Subroutine
