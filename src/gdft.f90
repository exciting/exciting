!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gdft (ik, delta)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Real (8), Intent (Out) :: delta (nstsv, nstsv)
! local variables
      Integer :: ist, jst, is, ia, ias
      Integer :: ir, irc, itp
      Real (8) :: sum, t1, t2
! automatic arrays
      Real (8) :: fr (nrcmtmax), gr (nrcmtmax), cf (3, nrcmtmax)
      Real (8) :: rflm (lmmaxvr), rftp (lmmaxvr)
! allocatable arrays
      Real (8), Allocatable :: rfmt (:, :, :)
      Real (8), Allocatable :: rfir (:)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: wfmt (:, :, :, :, :)
      Complex (8), Allocatable :: wfir (:, :, :)
      Allocate (rfmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (rfir(ngrtot))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (wfmt(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir(ngrtot, nspinor, nstsv))
! get the eigenvectors from file for input k-point
      Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
      Call getevecsv (vkl(:, ik), evecsv)
! find the matching coefficients
      Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), &
     & sfacgk(:, :, 1, ik), apwalm)
! calculate the wavefunctions for all second-variational states
      Call genwfsv (.False., ngk(1, ik), igkig(:, 1, ik), evalsv(:, &
     & ik), apwalm, evecfv, evecsv, wfmt, wfir)
      Do ist = 1, nstsv
         delta (ist, ist) = 0.d0
         Do jst = ist + 1, nstsv
! muffin-tin part
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  irc = 0
                  Do ir = 1, nrmt (is), input%groundstate%lradstep
                     irc = irc + 1
                     rflm (:) = 2.d0 * (exmt(:, ir, ias)+ecmt(:, ir, &
                    & ias)) - vxcmt (:, ir, ias)
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                    & lmmaxvr, rflm, 1, 0.d0, rftp, 1)
                     Do itp = 1, lmmaxvr
                        t1 = dble (wfmt(itp, irc, ias, 1, ist)) ** 2 + &
                       & aimag (wfmt(itp, irc, ias, 1, ist)) ** 2
                        t2 = dble (wfmt(itp, irc, ias, 1, jst)) ** 2 + &
                       & aimag (wfmt(itp, irc, ias, 1, jst)) ** 2
                        If (associated(input%groundstate%spin)) Then
                           t1 = t1 + dble (wfmt(itp, irc, ias, 2, ist)) &
                          & ** 2 + aimag (wfmt(itp, irc, ias, 2, ist)) &
                          & ** 2
                           t2 = t2 + dble (wfmt(itp, irc, ias, 2, jst)) &
                          & ** 2 + aimag (wfmt(itp, irc, ias, 2, jst)) &
                          & ** 2
                        End If
                        rftp (itp) = rftp (itp) * (t1-t2)
                     End Do
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                    & lmmaxvr, rftp, 1, 0.d0, rfmt(1, irc, ias), 1)
                  End Do
               End Do
            End Do
! interstitial part
            Do ir = 1, ngrtot
               t1 = dble (wfir(ir, 1, ist)) ** 2 + aimag (wfir(ir, 1, &
              & ist)) ** 2
               t2 = dble (wfir(ir, 1, jst)) ** 2 + aimag (wfir(ir, 1, &
              & jst)) ** 2
               If (associated(input%groundstate%spin)) Then
                  t1 = t1 + dble (wfir(ir, 2, ist)) ** 2 + aimag &
                 & (wfir(ir, 2, ist)) ** 2
                  t2 = t2 + dble (wfir(ir, 2, jst)) ** 2 + aimag &
                 & (wfir(ir, 2, jst)) ** 2
               End If
               rfir (ir) = (2.d0*(exir(ir)+ecir(ir))-vxcir(ir)) * &
              & (t1-t2)
            End Do
! integrate function
            sum = 0.d0
            Do ir = 1, ngrtot
               sum = sum + rfir (ir) * cfunir (ir)
            End Do
            sum = sum * omega / dble (ngrtot)
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  Do irc = 1, nrcmt (is)
                     fr (irc) = rfmt (1, irc, ias) * rcmt (irc, is) ** &
                    & 2
                  End Do
                  Call fderiv (-1, nrcmt(is), rcmt(1, is), fr, gr, cf)
                  sum = sum + fourpi * y00 * gr (nrcmt(is))
               End Do
            End Do
            delta (ist, jst) = sum
            delta (jst, ist) = - sum
         End Do
      End Do
      Deallocate (rfmt, rfir, wfmt, wfir)
      Deallocate (apwalm, evecfv, evecsv)
      Return
End Subroutine
