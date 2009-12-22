!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine exxengyk (ikp, evv, ecv)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Real (8), Intent (Inout) :: evv
      Real (8), Intent (Inout) :: ecv
! local variables
      Integer :: ngknr, ik, ist, jst
      Integer :: is, ia, ias, nrc, m, lmax
      Integer :: iv (3), iq, ig, igq0
      Real (8) :: cfq, v (3), t1
      Complex (8) zrho0, zt1
! automatic arrays
      Real (8) :: zn (nspecies)
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
      Real (8), Allocatable :: evalsvp (:)
      Real (8), Allocatable :: evalsvnr (:)
      Complex (8), Allocatable :: sfacgknr (:, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :, :, :)
      Complex (8), Allocatable :: wfir1 (:, :, :)
      Complex (8), Allocatable :: wfir2 (:, :, :)
      Complex (8), Allocatable :: wfcr (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :)
      Complex (8), Allocatable :: zvclir (:)
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
      Allocate (evalsvp(nstsv))
      Allocate (evalsvnr(nstsv))
      Allocate (sfacgknr(ngkmax, natmtot))
      Allocate (ylmgq(lmmaxvr, ngvec))
      Allocate (sfacgq(ngvec, natmtot))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir1(ngrtot, nspinor, nstsv))
      Allocate (wfir2(ngrtot, nspinor, nstsv))
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zvclir(ngrtot))
      Allocate (wfcr(lmmaxvr, nrcmtmax, 2))
      Allocate (zfmt(lmmaxvr, nrcmtmax))
! coefficient for long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the nuclear charges to zero
      zn (:) = 0.d0
! get the eigenvalues/vectors from file for input k-point
      Call getevalsv (vkl(:, ikp), evalsvp)
      Call getevecfv (vkl(:, ikp), vgkl(:, :, :, ikp), evecfv)
      Call getevecsv (vkl(:, ikp), evecsv)
! find the matching coefficients
      Call match (ngk(1, ikp), gkc(:, 1, ikp), tpgkc(:, :, 1, ikp), &
     & sfacgk(:, :, 1, ikp), apwalm)
! calculate the wavefunctions for occupied states for the input k-point
      Call genwfsv (.True., ngk(1, ikp), igkig(:, 1, ikp), evalsvp, &
     & apwalm, evecfv, evecsv, wfmt1, wfir1)
! start loop over non-reduced k-point set
      Do ik = 1, nkptnr
! generate G+k vectors
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
! spherical harmonics for G+q-vectors
            Call genylm (input%groundstate%lmaxvr, tpgqc(:, ig), &
           & ylmgq(:, ig))
         End Do
! structure factor for G+q
         Call gensfacgp (ngvec, vgqc, ngvec, sfacgq)
! find the shortest G+q-vector
         Call findigp0 (ngvec, gqc, igq0)
! compute the required spherical Bessel functions
         lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
         Call genjlgpr (lmax, gqc, jlgqr)
! calculate the wavefunctions for occupied states
         Call genwfsv (.True., ngknr, igkignr, evalsvnr, apwalm, &
        & evecfv, evecsv, wfmt2, wfir2)
!--------------------------------------------!
!    valence-valence-valence contribution    !
!--------------------------------------------!
         Do jst = 1, nstsv
            If (evalsvnr(jst) .Lt. efermi) Then
               Do ist = 1, nstsv
                  If (evalsvp(ist) .Lt. efermi) Then
! calculate the complex overlap density
                     Call vnlrho (.True., wfmt2(:, :, :, :, jst), &
                    & wfmt1(:, :, :, :, ist), wfir2(:, :, jst), &
                    & wfir1(:, :, ist), zrhomt, zrhoir)
! calculate the Coulomb potential
                     Call zpotcoul (nrcmt, nrcmtmax, nrcmtmax, rcmt, &
                    & igq0, gqc, jlgqr, ylmgq, sfacgq, zn, zrhomt, &
                    & zrhoir, zvclmt, zvclir, zrho0)
                     zt1 = zfinp (.True., zrhomt, zvclmt, zrhoir, &
                    & zvclir)
                     t1 = cfq * wiq2 (iq) * &
                    & (dble(zrho0)**2+aimag(zrho0)**2)
!$OMP CRITICAL
                     evv = evv - 0.5d0 * occmax * wkpt (ikp) * &
                    & (wkptnr(ik)*dble(zt1)+t1)
!$OMP END CRITICAL
! end loop over ist
                  End If
               End Do
! end loop over jst
            End If
         End Do
! end loop over non-reduced k-point set
      End Do
!-----------------------------------------!
!    valence-core-valence contribution    !
!-----------------------------------------!
! begin loops over atoms and species
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do jst = 1, spnst (is)
               If (spcore(jst, is)) Then
                  Do m = - spk (jst, is), spk (jst, is) - 1
! pass m-1/2 to wavefcr
                     Call wavefcr (input%groundstate%lradstep, is, ia, &
                    & jst, m, nrcmtmax, wfcr)
                     Do ist = 1, nstsv
                        If (evalsvp(ist) .Lt. efermi) Then
! calculate the complex overlap density
                           Call vnlrhomt (.True., is, wfcr(:, :, 1), &
                          & wfmt1(:, :, ias, 1, ist), zrhomt(:, :, &
                          & ias))
                           If (associated(input%groundstate%spin)) Then
                              Call vnlrhomt (.True., is, wfcr(:, :, 2), &
                             & wfmt1(:, :, ias, 2, ist), zfmt)
                              zrhomt (:, 1:nrc, ias) = zrhomt (:, &
                             & 1:nrc, ias) + zfmt (:, 1:nrc)
                           End If
! calculate the Coulomb potential
                           Call zpotclmt (input%groundstate%ptnucl, &
                          & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                          & 0.d0, lmmaxvr, zrhomt(:, :, ias), zvclmt(:, &
                          & :, ias))
                           zt1 = zfmtinp (.True., &
                          & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                          & lmmaxvr, zrhomt(:, :, ias), zvclmt(:, :, &
                          & ias))
!$OMP CRITICAL
                           ecv = ecv - occmax * wkpt (ikp) * dble (zt1)
!$OMP END CRITICAL
! end loop over ist
                        End If
                     End Do
! end loop over m
                  End Do
! end loop over jst
               End If
            End Do
! end loops over atoms and species
         End Do
      End Do
      Deallocate (igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
      Deallocate (vgqc, tpgqc, gqc, jlgqr)
      Deallocate (evalsvp, evalsvnr, evecfv, evecsv)
      Deallocate (sfacgknr, apwalm, ylmgq, sfacgq)
      Deallocate (wfmt1, wfmt2, wfir1, wfir2, wfcr)
      Deallocate (zrhomt, zrhoir, zvclmt, zvclir, zfmt)
      Return
End Subroutine
