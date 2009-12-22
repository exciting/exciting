!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phdveff (iph, iq, veffmtp, veffirp, dveffmt, dveffir)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iph
      Integer, Intent (In) :: iq
      Real (8), Intent (In) :: veffmtp (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (In) :: veffirp (ngrtot)
      Complex (8), Intent (Out) :: dveffmt (lmmaxvr, nrcmtmax, &
     & natmtot0)
      Complex (8), Intent (Out) :: dveffir (ngrtot0)
! local variables
      Integer :: is, ia, ja, ias, jas
      Integer :: ir, irc, i1, i2, i3, i
      Real (8) :: v1 (3), v2 (3), v3 (3), t1, t2
      Complex (8) zt1, zt2
! automatic arrays
      Real (8) :: rflm (lmmaxvr)
      Complex (8) zflm (lmmaxvr)
! external functions
      Real (8) :: rfirvec
      External rfirvec
! prefactor
      zt1 = 1.d0 / (dble(nphcell)*input%phonons%deltaph)
! multiply by i for sin-like displacement
      If (iph .Eq. 1) zt1 = zt1 * zi
!------------------------------!
!     muffin-tin potential     !
!------------------------------!
      ias = 0
      jas = 0
      Do is = 1, nspecies
         ja = 0
         Do ia = 1, natoms0 (is)
            ias = ias + 1
            Do i = 1, nphcell
               ja = ja + 1
               jas = jas + 1
! important: the muffin-tin potential should have an *explicit* phase exp(iq.r)
               t1 = - dot_product (vqc(:, iq), atposc(:, ja, is))
               zt2 = zt1 * cmplx (Cos(t1), Sin(t1), 8)
! loop over radial points
               irc = 0
               Do ir = 1, nrmt (is), input%groundstate%lradstep
                  irc = irc + 1
! compute the difference between the perturbed and unperturbed potentials
                  rflm (:) = veffmt (:, ir, jas) - veffmtp (:, ir, jas)
! convert real potential to a complex spherical harmonic expansion
                  Call rtozflm (input%groundstate%lmaxvr, rflm, zflm)
! add to total
                  dveffmt (:, irc, ias) = dveffmt (:, irc, ias) + zt2 * &
                 & zflm (:)
! end loop over radial points
               End Do
            End Do
! end loop over atoms and species
         End Do
      End Do
!--------------------------------!
!     interstitial potential     !
!--------------------------------!
      ir = 0
      Do i3 = 0, ngrid0 (3) - 1
         v1 (3) = dble (i3) / dble (ngrid0(3))
         Do i2 = 0, ngrid0 (2) - 1
            v1 (2) = dble (i2) / dble (ngrid0(2))
            Do i1 = 0, ngrid0 (1) - 1
               v1 (1) = dble (i1) / dble (ngrid0(1))
               ir = ir + 1
               Call r3mv (avec0, v1, v2)
               Do i = 1, nphcell
                  v3 (:) = v2 (:) + vphcell (:, i)
                  t1 = - dot_product (vqc(:, iq), v3(:))
                  zt2 = zt1 * cmplx (Cos(t1), Sin(t1), 8)
                  t1 = rfirvec (ngrid, ainv, v3, veffir)
                  t2 = rfirvec (ngrid, ainv, v3, veffirp)
                  dveffir (ir) = dveffir (ir) + zt2 * (t1-t2)
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
