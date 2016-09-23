!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gwf2cr (gwf2mt)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (Inout) :: gwf2mt (lmmaxvr, nrmtmax, natmtot)
! local variables
      Integer :: is, ia, ias, ist, i
      Integer :: l, m, lm, ir, itp
      Real (8) :: t1
! automatic arrays
      Complex (8) zftp (lmmaxvr)
! allocatable arrays
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: gzfmt (:, :, :)
      Allocate (zfmt(lmmaxvr, nrmtmax))
      Allocate (gzfmt(lmmaxvr, nrmtmax, 3))
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, spnst (is)
               If (spcore(ist, is) .And. (spk(ist, is) .Eq. spl(ist, &
              & is)+1)) Then
                  l = spl (ist, is)
                  Do m = - l, l
                     lm = idxlm (l, m)
                     zfmt (:, 1:nrmt(is)) = 0.d0
                     Do ir = 1, nrmt (is)
                        zfmt (lm, ir) = rwfcr (ir, 1, ist, ias) / spr &
                       & (ir, is)
                     End Do
                     Call gradzfmt (input%groundstate%lmaxvr, nrmt(is), &
                    & spr(:, is), lmmaxvr, nrmtmax, zfmt, gzfmt)
                     Do i = 1, 3
                        Do ir = 1, nrmt (is)
                           Call zgemv ('N', lmmaxvr, lmmaxvr, zone, &
                          & zbshtvr, lmmaxvr, gzfmt(:, ir, i), 1, &
                          & zzero, zftp, 1)
                           Do itp = 1, lmmaxvr
                              t1 = dble (zftp(itp)) ** 2 + aimag &
                             & (zftp(itp)) ** 2
! factor of 2 from spin
                              gwf2mt (itp, ir, ias) = gwf2mt (itp, ir, &
                             & ias) + 2.d0 * t1
                           End Do
                        End Do
                     End Do
                  End Do
               End If
            End Do
! end loops over atoms and species
         End Do
      End Do
      Deallocate (zfmt, gzfmt)
      Return
End Subroutine
