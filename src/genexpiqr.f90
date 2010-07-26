
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine genexpiqr (ik, emat)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (Out) :: emat (nstsv, nstsv)
! local variables
      Integer :: l1, l2, l3, m1, m2, m3
      Integer :: lm1, lm2, lm3, ist, jst
      Integer :: is, ia, i1, i2, i3, iv (3)
      Integer :: ngkq, igk, ifg, ir, irc
      Integer :: i, j, k, l, ispn
      Real (8) :: vecqc (3), qc, tp (2)
      Real (8) :: vkql (3), vkqc (3), x, t1
      Real (8) :: v1 (3), v2 (3), v3 (3)
      Complex (8) zsum, zt1, zt2, zt3
! automatic arrays
      Complex (8) ylm (lmmaxvr), zl (0:input%groundstate%lmaxvr)
      Complex (8) zflm (lmmaxvr)
! allocatable arrays
      Integer, Allocatable :: igkqig (:)
      Real (8), Allocatable :: gnt (:, :, :)
      Real (8), Allocatable :: jlqr (:, :)
      Real (8), Allocatable :: vgkql (:, :)
      Real (8), Allocatable :: vgkqc (:, :)
      Real (8), Allocatable :: gkqc (:)
      Real (8), Allocatable :: tpgkqc (:, :)
      Complex (8), Allocatable :: sfacgkq (:, :)
      Complex (8), Allocatable :: apwalm1 (:, :, :, :)
      Complex (8), Allocatable :: apwalm2 (:, :, :, :)
      Complex (8), Allocatable :: evecfv1 (:, :)
      Complex (8), Allocatable :: evecfv2 (:, :)
      Complex (8), Allocatable :: evecsv1 (:, :)
      Complex (8), Allocatable :: evecsv2 (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Complex (8), Allocatable :: wfmt3 (:, :)
      Complex (8), Allocatable :: wfir (:)
      Complex (8), Allocatable :: zfir1 (:)
      Complex (8), Allocatable :: zfir2 (:)
      Complex (8), Allocatable :: em (:, :)
! external functions
      Real (8) :: gaunt
      Complex (8) zfmtinp, zdotc
      External gaunt, zfmtinp, zdotc
! check if q-vector is zero
      t1 = input%properties%elnes%vecql(1) ** 2 + &
     & input%properties%elnes%vecql(2) ** 2 + &
     & input%properties%elnes%vecql(3) ** 2
      If (t1 .Lt. input%structure%epslat) Then
         emat (:, :) = 0.d0
         Do i = 1, nstsv
            emat (i, i) = 1.d0
         End Do
         Return
      End If
! check q-vector is commensurate with k-point grid
      v1 (:) = dble (input%groundstate%ngridk(:)) * &
     & input%properties%elnes%vecql(:)
      v2 (:) = Abs (v1(:)-Nint(v1(:)))
      If ((v2(1) .Gt. input%structure%epslat) .Or. (v2(2) .Gt. &
     & input%structure%epslat) .Or. (v2(3) .Gt. &
     & input%structure%epslat)) Then
         Write (*,*)
         Write (*, '("Error(genexpiqr): q-vector incommensurate with k-&
        &point grid")')
         Write (*, '(" ngridk : ", 3I6)') input%groundstate%ngridk
         Write (*, '(" vecql : ", 3G18.10)') &
        & input%properties%elnes%vecql
         Write (*,*)
         Stop
      End If
! allocate local arrays
      Allocate (igkqig(ngkmax))
      Allocate (gnt(lmmaxvr, lmmaxvr, lmmaxvr))
      Allocate (jlqr(0:input%groundstate%lmaxvr, nrcmtmax))
      Allocate (vgkql(3, ngkmax))
      Allocate (vgkqc(3, ngkmax))
      Allocate (gkqc(ngkmax))
      Allocate (tpgkqc(2, ngkmax))
      Allocate (sfacgkq(ngkmax, natmtot))
      Allocate (apwalm1(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (apwalm2(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv1(nmatmax, nstfv))
      Allocate (evecfv2(nmatmax, nstfv))
      If (input%groundstate%tevecsv) Then
         Allocate (evecsv1(nstsv, nstsv))
         Allocate (evecsv2(nstsv, nstsv))
      End If
      Allocate (wfmt1(lmmaxvr, nrcmtmax))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, nstfv))
      Allocate (wfmt3(lmmaxvr, nrcmtmax))
      Allocate (wfir(ngkmax))
      Allocate (zfir1(ngrtot), zfir2(ngrtot))
      Allocate (em(nstfv, nstfv))
! compute the Gaunt coefficients
      Do l1 = 0, input%groundstate%lmaxvr
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do l2 = 0, input%groundstate%lmaxvr
               Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  Do l3 = 0, input%groundstate%lmaxvr
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        gnt (lm1, lm2, lm3) = gaunt (l1, l2, l3, m1, &
                       & m2, m3)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
! q-vector in Cartesian coordinates
      Call r3mv (bvec, input%properties%elnes%vecql, vecqc)
! length and spherical coordinates of q-vector
      Call sphcrd (vecqc, qc, tp)
! generate the conjugate spherical harmonics of the q-vector
      Call genylm (input%groundstate%lmaxvr, tp, ylm)
      ylm (:) = conjg (ylm(:))
! get the eigenvector for k-point k
      Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv1)
! find the matching coefficients for k-point k
      Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), &
     & sfacgk(:, :, 1, ik), apwalm1)
! k+q-vector in lattice coordinates
      vkql (:) = vkl (:, ik) + input%properties%elnes%vecql(:)
! map vector components to [0,1) interval
      Call r3frac (input%structure%epslat, vkql, iv)
! k+q-vector in Cartesian coordinates
      Call r3mv (bvec, vkql, vkqc)
! generate the G+k+q-vectors
      Call gengpvec (vkql, vkqc, ngkq, igkqig, vgkql, vgkqc, gkqc, &
     & tpgkqc)
! generate the structure factors
      Call gensfacgp (ngkq, vgkqc, ngkmax, sfacgkq)
! find the matching coefficients for k-point k+q
      Call match (ngkq, gkqc, tpgkqc, sfacgkq, apwalm2)
! get the eigenvector for k-point k+q
      Call getevecfv (vkql, vgkql, evecfv2)
! set the first-variational matrix element array to zero
      em (:, :) = 0.d0
!------------------------------------!
!     muffin-tin matrix elements     !
!------------------------------------!
      Do is = 1, nspecies
! compute the spherical Bessel functions
         Do irc = 1, nrcmt (is)
            x = qc * rcmt (irc, is)
            Call sbessel (input%groundstate%lmaxvr, x, jlqr(:, irc))
         End Do
         Do ia = 1, natoms (is)
            t1 = dot_product (vecqc(:), atposc(:, ia, is))
            zt1 = fourpi * cmplx (Cos(t1), Sin(t1), 8)
            Do l1 = 0, input%groundstate%lmaxvr
               zl (l1) = zt1 * zil (l1)
            End Do
            Do ist = 1, nstfv
! calculate the wavefunction for k-point k+q
               Call wavefmt (input%groundstate%lradstep, &
              & input%groundstate%lmaxvr, is, ia, ngkq, apwalm2, &
              & evecfv2(:, ist), lmmaxvr, wfmt2(:, :, ist))
            End Do
            Do jst = 1, nstfv
! calculate the wavefunction for k-point k
               Call wavefmt (input%groundstate%lradstep, &
              & input%groundstate%lmaxvr, is, ia, ngk(1, ik), apwalm1, &
              & evecfv1(:, jst), lmmaxvr, wfmt1)
! multiply wavefunction with exp(iq.r)
               Do irc = 1, nrcmt (is)
                  zflm (:) = 0.d0
                  Do l2 = 0, input%groundstate%lmaxvr
                     zt1 = zl (l2) * jlqr (l2, irc)
                     Do m2 = - l2, l2
                        lm2 = idxlm (l2, m2)
                        zt2 = zt1 * ylm (lm2)
                        Do lm3 = 1, lmmaxvr
                           zt3 = zt2 * wfmt1 (lm3, irc)
                           Do lm1 = 1, lmmaxvr
                              zflm (lm1) = zflm (lm1) + gnt (lm1, lm2, &
                             & lm3) * zt3
                           End Do
                        End Do
                     End Do
                  End Do
                  wfmt3 (:, irc) = zflm (:)
               End Do
               Do ist = 1, nstfv
                  em (ist, jst) = em (ist, jst) + zfmtinp (.True., &
                 & input%groundstate%lmaxvr, nrcmt(is), rcmt(:, is), &
                 & lmmaxvr, wfmt2(:, :, ist), wfmt3)
               End Do
            End Do
! end loops over atoms and species
         End Do
      End Do
!--------------------------------------!
!     interstitial matrix elements     !
!--------------------------------------!
! store q+k-k', where k' is the k+q-vector mapped to [0,1)
      v1 (:) = vecqc (:) + vkc (:, ik) - vkqc (:)
! compute exp(i(q+k-k').r) times by the characteristic function
      ir = 0
      Do i3 = 0, ngrid (3) - 1
         v2 (3) = dble (i3) / dble (ngrid(3))
         Do i2 = 0, ngrid (2) - 1
            v2 (2) = dble (i2) / dble (ngrid(2))
            Do i1 = 0, ngrid (1) - 1
               v2 (1) = dble (i1) / dble (ngrid(1))
               ir = ir + 1
               Call r3mv (input%structure%crystal%basevect, v2, v3)
               t1 = dot_product (v1(:), v3(:))
               zfir1 (ir) = cfunir (ir) * cmplx (Cos(t1), Sin(t1), 8)
            End Do
         End Do
      End Do
! compute interstitial wavefunctions for k-point k
      Do jst = 1, nstfv
         zfir2 (:) = 0.d0
         Do igk = 1, ngk (1, ik)
            ifg = igfft (igkig(igk, 1, ik))
            zfir2 (ifg) = evecfv1 (igk, jst)
         End Do
! Fourier transform wavefunction to real-space
         Call zfftifc (3, ngrid, 1, zfir2)
! multiply with the phase and characteristic function
         zfir2 (:) = zfir2 (:) * zfir1 (:)
! Fourier transform back to G-space
         Call zfftifc (3, ngrid,-1, zfir2)
! store in wfir
         Do igk = 1, ngkq
            ifg = igfft (igkqig(igk))
            wfir (igk) = zfir2 (ifg)
         End Do
! add to the first-variational matrix elements
         Do ist = 1, nstfv
            em (ist, jst) = em (ist, jst) + zdotc (ngkq, evecfv2(:, &
           & ist), 1, wfir, 1)
         End Do
      End Do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
      If (input%groundstate%tevecsv) Then
! get the second-variational eigenvectors
         Call getevecsv (vkl(:, ik), evecsv1)
         Call getevecsv (vkql, evecsv2)
         Do i = 1, nstsv
            Do j = 1, nstsv
               zsum = 0.d0
               k = 0
               Do ispn = 1, nspinor
                  Do ist = 1, nstfv
                     k = k + 1
                     l = (ispn-1) * nstfv
                     Do jst = 1, nstfv
                        l = l + 1
                        zsum = zsum + em (ist, jst) * conjg (evecsv2(k, &
                       & i)) * evecsv1 (l, j)
                     End Do
                  End Do
               End Do
               emat (i, j) = zsum
            End Do
         End Do
      Else
         emat (:, :) = em (:, :)
      End If
      Deallocate (igkqig, gnt, jlqr, vgkql, vgkqc, gkqc, tpgkqc)
      Deallocate (sfacgkq, apwalm1, apwalm2, evecfv1, evecfv2)
      If (input%groundstate%tevecsv) deallocate (evecsv1, evecsv2)
      Deallocate (wfmt1, wfmt2, wfmt3, wfir, zfir1, zfir2, em)
      Return
End Subroutine
