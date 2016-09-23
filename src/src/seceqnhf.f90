!
!
!
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine seceqnhf (ikp, evecsvp)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (Inout) :: evecsvp (nstsv, nstsv)
! local variables
      Integer :: is, ia, ias, ir, irc
      Integer :: ngknr, ik, jk, ist1, ist2, ist3
      Integer :: iq, ig, iv (3), igq0
      Integer :: lmax, lwork, info
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
      Real (8), Allocatable :: rwork (:)
      Real (8), Allocatable :: rfmt (:, :, :)
      Complex (8), Allocatable :: h (:, :)
      Complex (8), Allocatable :: vmat (:, :)
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
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: work (:)
! external functions
      Complex (8) zfinp
      External zfinp
!$OMP CRITICAL
      Write (*, '("Info(seceqnhf): ", I6, " of ", I6, " k-points")') &
     & ikp, nkpt
!$OMP END CRITICAL
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
      Allocate (rwork(3*nstsv))
      Allocate (h(nstsv, nstsv))
      Allocate (vmat(nstsv, nstsv))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (sfacgknr(ngkmax, natmtot))
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
      Allocate (zfmt(lmmaxvr, nrcmtmax))
      lwork = 2 * nstsv
      Allocate (work(lwork))
      Allocate (rfmt(lmmaxvr, nrcmtmax, natmtot))
! coefficient of long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the point charges to zero
      zn (:) = 0.d0
! get the eigenvalues/vectors from file for input k-point
      Call getevalsv (vkl(:, ikp), evalsvp)
      Call getevecfv (vkl(:, ikp), vgkl(:, :, :, ikp), evecfv)
! find the matching coefficients
      Call match (ngk(1, ikp), gkc(:, 1, ikp), tpgkc(:, :, 1, ikp), &
     & sfacgk(:, :, 1, ikp), apwalm)
! calculate the wavefunctions for all states for the input k-point
      Call genwfsv (.False., ngk(1, ikp), igkig(:, 1, ikp), evalsvp, &
     & apwalm, evecfv, evecsvp, wfmt1, wfir1)
! compute the new kinetic matrix elements
      Call zgemm ('N', 'N', nstsv, nstsv, nstsv, zone, kinmatc(:, :, &
     & ikp), nstsv, evecsvp, nstsv, zzero, vmat, nstsv)
      Call zgemm ('C', 'N', nstsv, nstsv, nstsv, zone, evecsvp, nstsv, &
     & vmat, nstsv, zzero, h, nstsv)
! convert muffin-tin Coulomb potential to spherical coordinates
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, vclmt(:, ir, ias), 1, 0.d0, rfmt(:, irc, ias), &
              & 1)
            End Do
         End Do
      End Do
! compute the Coulomb matrix elements and add
      Call genvmatk (rfmt, vclir, wfmt1, wfir1, vmat)
      h (:, :) = h (:, :) + vmat (:, :)
! zero the non-local matrix elements for passed k-point
      vmat (:, :) = 0.d0
! start loop over non-reduced k-point set
      Do ik = 1, nkptnr
! find the equivalent reduced k-point
         iv (:) = ivknr (:, ik)
         jk = ikmap (iv(1), iv(2), iv(3))
! generate the G+k vectors
         Call gengpvec (vklnr(:, ik), vkcnr(:, ik), ngknr, igkignr, &
        & vgklnr, vgkcnr, gkcnr, tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-point
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
! calculate the wavefunctions for all states
         Call genwfsv (.False., ngknr, igkignr, evalsvnr, apwalm, &
        & evecfv, evecsv, wfmt2, wfir2)
         Do ist3 = 1, nstsv
            If (occsv(ist3, jk) .Gt. input%groundstate%epsocc) Then
               Do ist2 = 1, nstsv
! calculate the complex overlap density
                  Call vnlrho (.True., wfmt2(:, :, :, :, ist3), &
                 & wfmt1(:, :, :, :, ist2), wfir2(:, :, ist3), wfir1(:, &
                 & :, ist2), zrhomt, zrhoir)
! calculate the Coulomb potential
                  Call zpotcoul (nrcmt, nrcmtmax, nrcmtmax, rcmt, igq0, &
                 & gqc, jlgqr, ylmgq, sfacgq, zn, zrhomt, zrhoir, &
                 & zvclmt, zvclir, zrho02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
                  Do ist1 = 1, ist2
! calculate the complex overlap density
                     Call vnlrho (.True., wfmt2(:, :, :, :, ist3), &
                    & wfmt1(:, :, :, :, ist1), wfir2(:, :, ist3), &
                    & wfir1(:, :, ist1), zrhomt, zrhoir)
                     zt1 = zfinp (.True., zrhomt, zvclmt, zrhoir, &
                    & zvclir)
! compute the density coefficient of the smallest G+q-vector
                     Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, igq0), &
                    & sfacgq0, zrhomt, zrhoir, zrho01)
                     zt2 = cfq * wiq2 (iq) * (conjg(zrho01)*zrho02)
                     t1 = occsv (ist3, jk) / occmax
                     vmat (ist1, ist2) = vmat (ist1, ist2) - t1 * &
                    & (wkptnr(ik)*zt1+zt2)
                  End Do
               End Do
            End If
         End Do
! end loop over non-reduced k-point set
      End Do
! add the non-local matrix elements to Hamiltonian
      h (:, :) = h (:, :) + vmat (:, :)
! diagonalise the Hartree-Fock Hamiltonian (eigenvalues in global array)
      Call zheev ('V', 'U', nstsv, h, nstsv, evalsv(:, ikp), work, &
     & lwork, rwork, info)
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(seceqnhf): diagonalisation of the Hartree-F&
        &ock Hamiltonian failed")')
         Write (*, '(" for k-point ", I8)') ikp
         Write (*, '(" ZHEEV returned INFO = ", I8)') info
         Write (*,*)
         Stop
      End If
! apply unitary transformation to second-variational states
      evecsv (:, :) = evecsvp (:, :)
      Call zgemm ('N', 'N', nstsv, nstsv, nstsv, zone, evecsv, nstsv, &
     & h, nstsv, zzero, evecsvp, nstsv)
      Deallocate (igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
      Deallocate (evalsvp, evalsvnr, evecfv, evecsv, rwork)
      Deallocate (h, vmat, apwalm, sfacgknr, ylmgq, sfacgq)
      Deallocate (wfmt1, wfmt2, wfir1, wfir2)
      Deallocate (zrhomt, zrhoir, zvclmt, zvclir, zfmt, work, rfmt)
      Return
End Subroutine
