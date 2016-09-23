!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dynrtoq (vpl, dynr, dynp)
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: vpl (3)
      Complex (8), Intent (In) :: dynr (3*natmtot, 3*natmtot, &
     & ngridq(1)*ngridq(2)*ngridq(3))
      Complex (8), Intent (Out) :: dynp (3*natmtot, 3*natmtot)
! local variables
      Integer :: i1, i2, i3, ir, i, j, is, js, ia, ja, ii,jj
      Real (8) :: vrl(3), vrsl(3)
      Complex (8) :: zwght
      Complex (8) :: dynr_ij(3,3), dynp_ij(3,3)
!
      dynp (:, :) = 0.d0
! loop over atoms and species
      i = 1
      Do is = 1,nspecies
         Do ia = 1,natoms(is)
            j = 1
            Do js = 1,nspecies
               Do ja = 1,natoms(js)
! loop over R-vectors
                  ir = 0
                  dynp_ij(:,:) = (0.d0, 0.d0)
                  Do i3 = ngridq (3) / 2 - ngridq (3) + 1, ngridq (3) / 2
                     Do i2 = ngridq (2) / 2 - ngridq (2) + 1, ngridq (2) / 2
                        Do i1 = ngridq (1) / 2 - ngridq (1) + 1, ngridq (1) / 2
                           ir = ir + 1
                           dynr_ij(:,:) = dynr(i:i+2,j:j+2,ir)
                           vrl(1) = dble(i1)
                           vrl(2) = dble(i2)
                           vrl(3) = dble(i3)
                           vrsl(:) = vrl(:) + input%structure%speciesarray(js)%species%atomarray(ja)%atom%coord(:) &
                                           & - input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                           call ws_weight(vrl, vrsl, vpl, zwght)
                           dynp_ij(:,:) = dynp_ij(:,:) + zwght*dynr_ij(:,:)
                        End Do
                     End Do
                  End Do
                  dynp(i:i+2,j:j+2) = dynp_ij(:,:)
                  j = j + 3
               End Do
            End Do
            i = i + 3
         End Do
      End Do
! symmetrise the dynamical matrix
      Call dynsym (vpl, dynp)
      Return
End Subroutine
