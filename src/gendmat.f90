!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gendmat (tspndg, tlmdg, lmin, lmax, is, ia, ngp, apwalm, &
& evecfv, evecsv, ld, dmat)
      Use modinput
      Use modmain
      Implicit None
! arguments
      Logical, Intent (In) :: tspndg
      Logical, Intent (In) :: tlmdg
      Integer, Intent (In) :: lmin
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp (nspnfv)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot, nspnfv)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: dmat (ld, ld, nspinor, nspinor, &
     & nstsv)
! local variables
      Integer :: ispn, jspn, lmmax
      Integer :: l, m1, m2, lm1, lm2
      Integer :: i, j, n, ist, irc
      Real (8) :: t1, t2
      Complex (8) zt1
! automatic arrays
      Real (8) :: fr1 (nrcmtmax), fr2 (nrcmtmax)
      Real (8) :: gr (nrcmtmax), cf (3, nrcmtmax)
! allocatable arrays
      Logical, Allocatable :: done (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      If (lmin .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(gendmat): lmin < 0 : ", I8)') lmin
         Write (*,*)
         Stop
      End If
      If (lmax .Gt. input%groundstate%lmaxapw) Then
         Write (*,*)
         Write (*, '("Error(gendmat): lmax > lmaxapw : ", 2I8)') lmax, &
        & input%groundstate%lmaxapw
         Write (*,*)
         Stop
      End If
      lmmax = (lmax+1) ** 2
! allocate local arrays
      Allocate (done(nstfv, nspnfv))
      Allocate (wfmt1(lmmax, nrcmtmax, nstfv, nspnfv))
      Allocate (wfmt2(lmmax, nrcmtmax, nspinor))
! zero the density matrix
      dmat (:, :, :, :, :) = 0.d0
      n = lmmax * nrcmt (is)
      done (:, :) = .False.
! begin loop over second-variational states
      Do j = 1, nstsv
         If (input%groundstate%tevecsv) Then
! generate spinor wavefunction from second-variational eigenvectors
            wfmt2 (:, :, :) = 0.d0
            i = 0
            Do ispn = 1, nspinor
               If (isspinspiral()) Then
                  jspn = ispn
               Else
                  jspn = 1
               End If
               Do ist = 1, nstfv
                  i = i + 1
                  zt1 = evecsv (i, j)
                  If (Abs(dble(zt1))+Abs(aimag(zt1)) .Gt. &
                 & input%groundstate%epsocc) Then
                     If ( .Not. done(ist, jspn)) Then
                        Call wavefmt (input%groundstate%lradstep, lmax, &
                       & is, ia, ngp(jspn), apwalm(:, :, :, :, jspn), &
                       & evecfv(:, ist, jspn), lmmax, wfmt1(:, :, ist, &
                       & jspn))
                        done (ist, jspn) = .True.
                     End If
! add to spinor wavefunction
                     Call zaxpy (n, zt1, wfmt1(:, :, ist, jspn), 1, &
                    & wfmt2(:, :, ispn), 1)
                  End If
               End Do
            End Do
         Else
! spin-unpolarised wavefunction
            Call wavefmt (input%groundstate%lradstep, lmax, is, ia, &
           & ngp, apwalm, evecfv(:, j, 1), lmmax, wfmt2)
         End If
         Do ispn = 1, nspinor
            Do jspn = 1, nspinor
               If (tspndg .And. (ispn .Ne. jspn)) Go To 20
               Do l = lmin, lmax
                  Do m1 = - l, l
                     lm1 = idxlm (l, m1)
                     Do m2 = - l, l
                        lm2 = idxlm (l, m2)
                        If (tlmdg .And. (lm1 .Ne. lm2)) Go To 10
                        Do irc = 1, nrcmt (is)
                           zt1 = wfmt2 (lm1, irc, ispn) * conjg &
                          & (wfmt2(lm2, irc, jspn))
                           t1 = rcmt (irc, is) ** 2
                           fr1 (irc) = dble (zt1) * t1
                           fr2 (irc) = aimag (zt1) * t1
                        End Do
                        Call fderiv (-1, nrcmt(is), rcmt(:, is), fr1, &
                       & gr, cf)
                        t1 = gr (nrcmt(is))
                        Call fderiv (-1, nrcmt(is), rcmt(:, is), fr2, &
                       & gr, cf)
                        t2 = gr (nrcmt(is))
                        dmat (lm1, lm2, ispn, jspn, j) = cmplx (t1, t2, &
                       & 8)
10                      Continue
                     End Do
                  End Do
               End Do
20             Continue
            End Do
         End Do
! end loop over second-variational states
      End Do
      Deallocate (done, wfmt1, wfmt2)
      Return
End Subroutine
