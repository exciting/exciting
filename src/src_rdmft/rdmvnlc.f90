!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmvnlc (ikp, vnl)
! calculate non-local matrix elements for minimisation w.r.t. evecsv
      Use modinput
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (Out) :: vnl (nstsv, nstsv, nstsv, nkptnr)
! local variables
      Integer :: ngknr, ik, ist1, ist2, ist3
      Integer :: lmax, ig, iq, igq0, iv (3)
      Real (8) :: cfq, v (3), t1
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
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :)
      Complex (8), Allocatable :: zvclir (:)
! external functions
      Complex (8) zfinp
      External zfinp
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
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zvclir(ngrtot))
! factor for long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the point charges to zero
      zn (:) = 0.d0
! get the eigenvectors and values from file
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
! structure factors for G+q
         Call gensfacgp (ngvec, vgqc, ngvec, sfacgq)
! find the shortest G+q-vector
         Call findigp0 (ngvec, gqc, igq0)
         sfacgq0 (:) = sfacgq (igq0, :)
! compute the required spherical Bessel functions
         lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
         Call genjlgpr (lmax, gqc, jlgqr)
         Call genjlgq0r (gqc(igq0), jlgq0r)
! calculate the wavefunctions for all states for non-reduced k-point ik
         Call genwfsv (.False., ngknr, igkignr, evalsvnr, apwalm, &
        & evecfv, evecsv, wfmt2, wfir2)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
         Do ist1 = 1, nstsv
            Do ist2 = 1, nstsv
! calculate the complex overlap density
               Call vnlrho (.True., wfmt2(:, :, :, :, ist2), wfmt1(:, &
              & :, :, :, ist1), wfir2(:, :, ist2), wfir1(:, :, ist1), &
              & zrhomt, zrhoir)
! compute the potential and G=0 coefficient of the density
               Call zpotcoul (nrcmt, nrcmtmax, nrcmtmax, rcmt, igq0, &
              & gqc, jlgqr, ylmgq, sfacgq, zn, zrhomt, zrhoir, zvclmt, &
              & zvclir, zrho02)
               zt1 = zfinp (.True., zrhomt, zvclmt, zrhoir, zvclir)
               t1 = cfq * wiq2 (iq) * &
              & (dble(zrho02)**2+aimag(zrho02)**2)
               vnl (ist1, ist1, ist2, ik) = wkptnr (ik) * dble (zt1) + &
              & t1
               Do ist3 = 1, nstsv
                  If (ist1 .Gt. ist3) Then
! calculate the complex overlap density
                     Call vnlrho (.True., wfmt2(:, :, :, :, ist2), &
                    & wfmt1(:, :, :, :, ist3), wfir2(:, :, ist2), &
                    & wfir1(:, :, ist3), zrhomt, zrhoir)
                     zt1 = zfinp (.True., zrhomt, zvclmt, zrhoir, &
                    & zvclir)
! compute the density coefficient of the smallest G+q-vector
                     Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, igq0), &
                    & sfacgq0, zrhomt, zrhoir, zrho01)
                     zt2 = cfq * wiq2 (iq) * (conjg(zrho01)*zrho02)
                     vnl (ist3, ist1, ist2, ik) = wkptnr (ik) * zt1 + &
                    & zt2
! end loop over ist3
                  End If
               End Do
! end loop over ist2
            End Do
! end loop over ist1
         End Do
! calculate the lower diagonal
         Do ist1 = 1, nstsv
            Do ist3 = 1, nstsv
               If (ist1 .Lt. ist3) Then
                  vnl (ist3, ist1, :, ik) = conjg (vnl(ist1, ist3, :, &
                 & ik))
               End If
            End Do
         End Do
! end loop over non-reduced k-point set
      End Do
      Deallocate (igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
      Deallocate (evalsvp, evalsvnr)
      Deallocate (apwalm, evecfv, evecsv, sfacgknr, ylmgq, sfacgq)
      Deallocate (wfmt1, wfmt2, wfir1, wfir2)
      Deallocate (zrhomt, zrhoir, zvclmt, zvclir)
      Return
End Subroutine
!EOC
