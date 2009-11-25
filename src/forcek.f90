!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: forcek
! !INTERFACE:
!
!
Subroutine forcek (ik, ffacg)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Computes the {\bf k}-dependent contribution to the incomplete basis set
!   (IBS) force. See the calling routine {\tt force} for a full description.
!
! !REVISION HISTORY:
!   Created June 2006 (JKD)
!   Updated for spin-spiral case, May 2007 (Francesco Cricchio and JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Real (8), Intent (In) :: ffacg (ngvec, nspecies)
! local variables
      Integer :: np, ispn, jspn
      Integer :: is, ia, ias, ist, jst
      Integer :: i, j, k, l, iv (3), ig
      Real (8) :: sum, t1
      Complex (8) zt1, zt2
      Complex (8) v (1)
! allocatable arrays
      Integer, Allocatable :: ijg (:)
      Real (8), Allocatable :: dp (:)
      Real (8), Allocatable :: evalfv (:, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: h (:)
      Complex (8), Allocatable :: o (:)
      Complex (8), Allocatable :: dlh (:)
      Complex (8), Allocatable :: dlo (:)
      Complex (8), Allocatable :: vh (:)
      Complex (8), Allocatable :: vo (:)
      Complex (8), Allocatable :: ffv (:, :)
      Complex (8), Allocatable :: y (:)
! external functions
      Complex (8) zdotc
      External zdotc
      np = npmat (1, ik)
      If (isspinspiral()) np = Max (np, npmat(2, ik))
! allocate local arrays
      Allocate (ijg(np))
      Allocate (dp(np))
      Allocate (evalfv(nstfv, nspnfv))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (h(np), o(np))
      Allocate (dlh(np), dlo(np))
      Allocate (vh(nmatmax))
      Allocate (vo(nmatmax))
      Allocate (ffv(nstfv, nstfv))
      Allocate (y(nstfv))
! get the eigenvalues/vectors and occupancies from file
      Call getevalfv (vkl(:, ik), evalfv)
      Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
      Call getevecsv (vkl(:, ik), evecsv)
      Call getoccsv (vkl(:, ik), occsv(:, ik))
! begin loop over first-variational spin components
      Do ispn = 1, nspnfv
! find the matching coefficients
         Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, &
        & ik), sfacgk(:, :, ispn, ik), apwalm)
         Do j = 1, ngk (ispn, ik)
            k = ((j-1)*j) / 2
            Do i = 1, j
               k = k + 1
               iv (:) = ivg (:, igkig(i, ispn, ik)) - ivg (:, igkig(j, &
              & ispn, ik))
               iv (:) = modulo (iv(:)-intgv(:, 1), ngrid(:)) + intgv &
              & (:, 1)
               ijg (k) = ivgig (iv(1), iv(2), iv(3))
               dp (k) = 0.5d0 * dot_product (vgkc(:, i, ispn, ik), &
              & vgkc(:, j, ispn, ik))
            End Do
         End Do
! loop over species and atoms
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               h (1:npmat(ispn, ik)) = 0.d0
               Call hmlaa (.False., is, ia, ngk(ispn, ik), apwalm, v, &
              & h)
               Call hmlalo (.False., is, ia, ngk(ispn, ik), apwalm, v, &
              & h)
               o (1:npmat(ispn, ik)) = 0.d0
               Call olpaa (.False., is, ia, ngk(ispn, ik), apwalm, v, &
              & o)
               Call olpalo (.False., is, ia, ngk(ispn, ik), apwalm, v, &
              & o)
! loop over Cartesian directions
               Do l = 1, 3
! APW-APW contribution
                  Do j = 1, ngk (ispn, ik)
                     k = ((j-1)*j) / 2
                     Do i = 1, j
                        k = k + 1
                        ig = ijg (k)
                        t1 = vgc (l, ig)
                        zt1 = - ffacg (ig, is) * conjg (sfacg(ig, ias))
                        dlh (k) = (dp(k)*zt1+h(k)) * t1
                        dlo (k) = (zt1+o(k)) * t1
                     End Do
                  End Do
! APW-local-orbital contribution
                  Do j = ngk (ispn, ik) + 1, nmat (ispn, ik)
                     k = ((j-1)*j) / 2
                     Do i = 1, ngk (ispn, ik)
                        k = k + 1
                        t1 = vgkc (l, i, ispn, ik)
                        dlh (k) = h (k) * t1
                        dlo (k) = o (k) * t1
                     End Do
                     Do i = ngk (ispn, ik) + 1, j
                        k = k + 1
                        dlh (k) = 0.d0
                        dlo (k) = 0.d0
                     End Do
                  End Do
! multiply by i
                  Do k = 1, npmat (ispn, ik)
                     dlh (k) = cmplx (-aimag(dlh(k)), dble(dlh(k)), 8)
                     dlo (k) = cmplx (-aimag(dlo(k)), dble(dlo(k)), 8)
                  End Do
! compute the force matrix elements in the first-variational basis
                  Do jst = 1, nstfv
                     Call zhpmv ('U', nmat(ispn, ik), zone, dlh, &
                    & evecfv(:, jst, ispn), 1, zzero, vh, 1)
                     Call zhpmv ('U', nmat(ispn, ik), zone, dlo, &
                    & evecfv(:, jst, ispn), 1, zzero, vo, 1)
                     t1 = evalfv (jst, ispn)
                     Do ist = 1, nstfv
                        zt1 = zdotc (nmat(ispn, ik), evecfv(:, ist, &
                       & ispn), 1, vh, 1)
                        zt2 = zdotc (nmat(ispn, ik), evecfv(:, ist, &
                       & ispn), 1, vo, 1)
                        ffv (ist, jst) = zt1 - t1 * zt2
                     End Do
                  End Do
! compute the force using the second-variational coefficients if required
                  sum = 0.d0
                  If (input%groundstate%tevecsv) Then
                     If (isspinspiral()) Then
! spin-spiral case
                        Do j = 1, nstsv
                           t1 = occsv (j, ik)
                           i = (ispn-1) * nstfv + 1
                           Call zgemv ('N', nstfv, nstfv, zone, ffv, &
                          & nstfv, evecsv(i, j), 1, zzero, y, 1)
                           zt1 = zdotc (nstfv, evecsv(i, j), 1, y, 1)
                           sum = sum + t1 * dble (zt1)
                        End Do
                     Else
! normal spin-polarised case
                        Do j = 1, nstsv
                           t1 = occsv (j, ik)
                           Do jspn = 1, nspinor
                              i = (jspn-1) * nstfv + 1
                              Call zgemv ('N', nstfv, nstfv, zone, ffv, &
                             & nstfv, evecsv(i, j), 1, zzero, y, 1)
                              zt1 = zdotc (nstfv, evecsv(i, j), 1, y, &
                             & 1)
                              sum = sum + t1 * dble (zt1)
                           End Do
                        End Do
                     End If
                  Else
! spin-unpolarised case
                     Do j = 1, nstsv
                        sum = sum + occsv (j, ik) * dble (ffv(j, j))
                     End Do
                  End If
!$OMP CRITICAL
                  forceibs (l, ias) = forceibs (l, ias) + wkpt (ik) * &
                 & sum
!$OMP END CRITICAL
! end loop over Cartesian components
               End Do
! end loop over atoms and species
            End Do
         End Do
! end loop over first-variational spins
      End Do
      Deallocate (ijg, dp, evalfv, apwalm, evecfv, evecsv)
      Deallocate (h, o, dlh, dlo, vh, vo, ffv, y)
      Return
End Subroutine
!EOC
