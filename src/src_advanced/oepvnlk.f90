!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine oepvnlk (ikp, vnlcv, vnlvv)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (Out) :: vnlcv (ncrmax, natmtot, nstsv)
      Complex (8), Intent (Out) :: vnlvv (nstsv, nstsv)
! local variables
      Integer :: ngknr, ik, ist1, ist2, ist3
      Integer :: is, ia, ias, ic, m1, m2, lmax
      Integer :: nrc, iq, ig, iv (3), igq0
      Real (8) :: v (3), cfq
      Complex (8) zrho01, zrho02, zt1, zt2
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Integer, Allocatable :: igkignr (:)
      Real (8), Allocatable :: vgklnr (:, :)
      Real (8), Allocatable :: vgkcnr (:, :)
      Real (8), Allocatable :: gkcnr (:)
      Real (8), Allocatable :: tpgkcnr (:, :)
      Real (8), Allocatable :: vgqc (:, :)
      Real (8), Allocatable :: tpgqc (:, :)
      Real (8), Allocatable :: gqc (:)
      Real (8), Allocatable :: jlgqr (:, :, :)
      Real (8), Allocatable :: jlgq0r (:, :, :)
      Real (8), Allocatable :: evalsvp (:)
      Real (8), Allocatable :: evalsvnr (:)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: sfacgknr (:, :)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :, :, :)
      Complex (8), Allocatable :: wfir1 (:, :, :)
      Complex (8), Allocatable :: wfir2 (:, :, :)
      Complex (8), Allocatable :: wfcr1 (:, :, :)
      Complex (8), Allocatable :: wfcr2 (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :)
      Complex (8), Allocatable :: zvclir (:)
      Complex (8), Allocatable :: zvcltp (:, :)
      Complex (8), Allocatable :: zfmt (:, :)
! external functions
      Complex (8) zfinp, zfmtinp
      External zfinp, zfmtinp
! allocate local arrays
      Allocate (igkignr(ngkmax))
      Allocate (vgklnr(3, ngkmax))
      Allocate (vgkcnr(3, ngkmax))
      Allocate (gkcnr(ngkmax))
      Allocate (tpgkcnr(2, ngkmax))
      Allocate (vgqc(3, ngvec))
      Allocate (tpgqc(2, ngvec))
      Allocate (gqc(ngvec))
      Allocate &
     & (jlgqr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, &
     & ngvec, nspecies))
      Allocate (jlgq0r(0:input%groundstate%lmaxvr, nrcmtmax, nspecies))
      Allocate (evalsvp(nstsv))
      Allocate (evalsvnr(nstsv))
      Allocate (sfacgknr(ngkmax, natmtot))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (ylmgq(lmmaxvr, ngvec))
      Allocate (sfacgq(ngvec, natmtot))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir1(ngrtot, nspinor, nstsv))
      Allocate (wfir2(ngrtot, nspinor, nstsv))
      Allocate (wfcr1(lmmaxvr, nrcmtmax, 2))
      Allocate (wfcr2(lmmaxvr, nrcmtmax, 2))
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zvclir(ngrtot))
      Allocate (zvcltp(lmmaxvr, nrcmtmax))
      Allocate (zfmt(lmmaxvr, nrcmtmax))
! factor for long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the nuclear charges to zero
      zn (:) = 0.d0
      vnlcv (:, :, :) = 0.d0
      vnlvv (:, :) = 0.d0
! get the eigenvalues/vectors from file for input k-point
      Call getevalsv (vkl(:, ikp), evalsvp)
      Call getevecfv (vkl(:, ikp), vgkl(:, :, :, ikp), evecfv)
      Call getevecsv (vkl(:, ikp), evecsv)
! find the matching coefficients
      Call match (ngk(1, ikp), gkc(:, 1, ikp), tpgkc(:, :, 1, ikp), &
     & sfacgk(:, :, 1, ikp), apwalm)
! calculate the wavefunctions for all states for the input k-point
      Call genwfsv (.False., ngk(1, ikp), igkig(:, 1, ikp), evalsvp, &
     & apwalm, evecfv, evecsv, wfmt1, wfir1)
! start loop over non-reduced k-point set
      Do ik = 1, nkptnr
! generate G+k-vectors
         Call gengpvec (vklnr(:, ik), vkcnr(:, ik), ngknr, igkignr, &
        & vgklnr, vgkcnr, gkcnr, tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-points
         Call getevalsv (vklnr(:, ik), evalsvnr)
         Call getevecfv (vklnr(:, ik), vgklnr, evecfv)
         Call getevecsv (vklnr(:, ik), evecsv)
! generate the structure factors
         Call gensfacgp (ngknr, vgkcnr, ngkmax, sfacgknr)
! find the matching coefficients
         Call match (ngknr, gkcnr, tpgkcnr, sfacgknr, apwalm)
! determine q-vector
         iv (:) = ivk (:, ikp) - ivknr (:, ik)
         iv (:) = modulo (iv(:), input%groundstate%ngridk(:))
         iq = iqmap (iv(1), iv(2), iv(3))
         v (:) = vkc (:, ikp) - vkcnr (:, ik)
         Do ig = 1, ngvec
! determine G+q vectors
            vgqc (:, ig) = vgc (:, ig) + v (:)
! G+q-vector length and (theta, phi) coordinates
            Call sphcrd (vgqc(:, ig), gqc(ig), tpgqc(:, ig))
! spherical harmonics for G+q-vector
            Call genylm (input%groundstate%lmaxvr, tpgqc(:, ig), &
           & ylmgq(:, ig))
         End Do
! structure factors for G+q
         Call gensfacgp (ngvec, vgqc, ngvec, sfacgq)
! find the shortest G+q-vector
         Call findigp0 (ngvec, gqc, igq0)
         sfacgq0 (:) = sfacgq (igq0, :)
! compute the required spherical Bessel functions
         lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
         Call genjlgpr (lmax, gqc, jlgqr)
         Call genjlgq0r (gqc(igq0), jlgq0r)
! calculate the wavefunctions for occupied states
         Call genwfsv (.True., ngknr, igkignr, evalsvnr, apwalm, &
        & evecfv, evecsv, wfmt2, wfir2)
         Do ist3 = 1, nstsv
            If (evalsvnr(ist3) .Lt. efermi) Then
               Do ist2 = 1, nstsv
                  If (evalsvp(ist2) .Gt. efermi) Then
! calculate the complex overlap density
                     Call vnlrho (.True., wfmt2(:, :, :, :, ist3), &
                    & wfmt1(:, :, :, :, ist2), wfir2(:, :, ist3), &
                    & wfir1(:, :, ist2), zrhomt, zrhoir)
! calculate the Coulomb potential
                     Call zpotcoul (nrcmt, nrcmtmax, nrcmtmax, rcmt, &
                    & igq0, gqc, jlgqr, ylmgq, sfacgq, zn, zrhomt, &
                    & zrhoir, zvclmt, zvclir, zrho02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
                     Do ist1 = 1, nstsv
                        If (evalsvp(ist1) .Lt. efermi) Then
! calculate the complex overlap density
                           Call vnlrho (.True., wfmt2(:, :, :, :, &
                          & ist3), wfmt1(:, :, :, :, ist1), wfir2(:, :, &
                          & ist3), wfir1(:, :, ist1), zrhomt, zrhoir)
                           zt1 = zfinp (.True., zrhomt, zvclmt, zrhoir, &
                          & zvclir)
! compute the density coefficient of the smallest G+q-vector
                           Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                          & igq0), sfacgq0, zrhomt, zrhoir, zrho01)
                           zt2 = cfq * wiq2 (iq) * &
                          & (conjg(zrho01)*zrho02)
                           vnlvv (ist1, ist2) = vnlvv (ist1, ist2) - &
                          & (wkptnr(ik)*zt1+zt2)
                        End If
                     End Do
!-------------------------------------------!
!     core-valence-valence contribution     !
!-------------------------------------------!
                     Do is = 1, nspecies
                        nrc = nrcmt (is)
                        Do ia = 1, natoms (is)
                           ias = idxas (ia, is)
! convert the muffin-tin potential to spherical coordinates
                           Call zgemm ('N', 'N', lmmaxvr, nrc, lmmaxvr, &
                          & zone, zbshtvr, lmmaxvr, zvclmt(:, :, ias), &
                          & lmmaxvr, zzero, zvcltp, lmmaxvr)
                           ic = 0
                           Do ist1 = 1, spnst (is)
                              If (spcore(ist1, is)) Then
                                 Do m1 = - spk (ist1, is), spk (ist1, &
                                & is) - 1
                                    ic = ic + 1
! pass m-1/2 to wavefcr
                                    Call wavefcr &
                                   & (input%groundstate%lradstep, is, &
                                   & ia, ist1, m1, nrcmtmax, wfcr1)
! calculate the complex overlap density
                                    Call vnlrhomt (.False., is, &
                                   & wfmt2(:, :, ias, 1, ist3), &
                                   & wfcr1(:, :, 1), zrhomt(:, :, ias))
                                    If (associated(input%groundstate%spin)) Then
                                       Call vnlrhomt (.False., is, &
                                      & wfmt2(:, :, ias, 2, ist3), &
                                      & wfcr1(:, :, 2), zfmt)
                                       zrhomt (:, 1:nrc, ias) = zrhomt &
                                      & (:, 1:nrc, ias) + zfmt (:, &
                                      & 1:nrc)
                                    End If
                                    zt1 = zfmtinp (.False., &
                                   & input%groundstate%lmaxvr, nrc, &
                                   & rcmt(:, is), lmmaxvr, zrhomt(:, :, &
                                   & ias), zvcltp)
                                    vnlcv (ic, ias, ist2) = vnlcv (ic, &
                                   & ias, ist2) - wkptnr (ik) * zt1
                                 End Do
! end loop over ist1
                              End If
                           End Do
! end loops over atoms and species
                        End Do
                     End Do
! end loop over ist2
                  End If
               End Do
! end loop over ist3
            End If
         End Do
! end loop over non-reduced k-point set
      End Do
! begin loops over atoms and species
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist3 = 1, spnst (is)
               If (spcore(ist3, is)) Then
                  Do m1 = - spk (ist3, is), spk (ist3, is) - 1
! pass m-1/2 to wavefcr
                     Call wavefcr (input%groundstate%lradstep, is, ia, &
                    & ist3, m1, nrcmtmax, wfcr1)
                     Do ist2 = 1, nstsv
                        If (evalsvp(ist2) .Gt. efermi) Then
! calculate the complex overlap density
                           Call vnlrhomt (.True., is, wfcr1(:, :, 1), &
                          & wfmt1(:, :, ias, 1, ist2), zrhomt(:, :, &
                          & ias))
                           If (associated(input%groundstate%spin)) Then
                              Call vnlrhomt (.True., is, wfcr1(:, :, &
                             & 2), wfmt1(:, :, ias, 2, ist2), zfmt)
                              zrhomt (:, 1:nrc, ias) = zrhomt (:, &
                             & 1:nrc, ias) + zfmt (:, 1:nrc)
                           End If
! calculate the Coulomb potential
                           Call zpotclmt (input%groundstate%ptnucl, &
                          & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                          & 0.d0, lmmaxvr, zrhomt(:, :, ias), zfmt)
! convert muffin-tin potential to spherical coordinates
                           Call zgemm ('N', 'N', lmmaxvr, nrc, lmmaxvr, &
                          & zone, zbshtvr, lmmaxvr, zfmt, lmmaxvr, &
                          & zzero, zvcltp, lmmaxvr)
!-------------------------------------------!
!     valence-core-valence contribution     !
!-------------------------------------------!
                           Do ist1 = 1, nstsv
                              If (evalsvp(ist1) .Lt. efermi) Then
! calculate the complex overlap density
                                 Call vnlrhomt (.False., is, wfcr1(:, &
                                & :, 1), wfmt1(:, :, ias, 1, ist1), &
                                & zrhomt(:, :, ias))
                                 If &
                                & (associated(input%groundstate%spin)) &
                                & Then
                                    Call vnlrhomt (.False., is, &
                                   & wfcr1(:, :, 2), wfmt1(:, :, ias, &
                                   & 2, ist1), zfmt)
                                    zrhomt (:, 1:nrc, ias) = zrhomt (:, &
                                   & 1:nrc, ias) + zfmt (:, 1:nrc)
                                 End If
                                 zt1 = zfmtinp (.False., &
                                & input%groundstate%lmaxvr, nrc, &
                                & rcmt(:, is), lmmaxvr, zrhomt(:, :, &
                                & ias), zvcltp)
                                 vnlvv (ist1, ist2) = vnlvv (ist1, &
                                & ist2) - zt1
                              End If
                           End Do
!----------------------------------------!
!     core-core-valence contribution     !
!----------------------------------------!
                           ic = 0
                           Do ist1 = 1, spnst (is)
                              If (spcore(ist1, is)) Then
                                 Do m2 = - spk (ist1, is), spk (ist1, &
                                & is) - 1
                                    ic = ic + 1
! pass m-1/2 to wavefcr
                                    Call wavefcr &
                                   & (input%groundstate%lradstep, is, &
                                   & ia, ist1, m2, nrcmtmax, wfcr2)
! calculate the complex overlap density
                                    Call vnlrhomt (.False., is, &
                                   & wfcr1(:, :, 1), wfcr2(:, :, 1), &
                                   & zrhomt(:, :, ias))
                                    Call vnlrhomt (.False., is, &
                                   & wfcr1(:, :, 2), wfcr2(:, :, 2), &
                                   & zfmt)
                                    zrhomt (:, 1:nrc, ias) = zrhomt (:, &
                                   & 1:nrc, ias) + zfmt (:, 1:nrc)
                                    zt1 = zfmtinp (.False., &
                                   & input%groundstate%lmaxvr, nrc, &
                                   & rcmt(:, is), lmmaxvr, zrhomt(:, :, &
                                   & ias), zvcltp)
                                    vnlcv (ic, ias, ist2) = vnlcv (ic, &
                                   & ias, ist2) - zt1
                                 End Do
! end loop over ist1
                              End If
                           End Do
! end loop over ist2
                        End If
                     End Do
! end loops over ist3 and m1
                  End Do
               End If
            End Do
! end loops over atoms and species
         End Do
      End Do
      Deallocate (igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
      Deallocate (evalsvp, evalsvnr, evecfv, evecsv)
      Deallocate (apwalm, sfacgknr, ylmgq, sfacgq)
      Deallocate (wfmt1, wfmt2, wfir1, wfir2, wfcr1, wfcr2)
      Deallocate (zrhomt, zrhoir, zvclmt, zvclir, zvcltp, zfmt)
      Return
End Subroutine
!EOC
