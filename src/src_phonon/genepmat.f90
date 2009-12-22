!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genepmat (iq, vpl, dveffmt, dveffir, epmat)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iq
      Real (8), Intent (In) :: vpl (3)
      Complex (8), Intent (In) :: dveffmt (lmmaxapw, nrcmtmax, natmtot, &
     & 3*natmtot)
      Complex (8), Intent (In) :: dveffir (ngrtot, 3*natmtot)
      Complex (8), Intent (Out) :: epmat (nstsv, nstsv, 3*natmtot)
! local variables
      Integer :: is, ia, ias
      Integer :: ngp, ngpq, igp, ifg
      Integer :: nrc, irc, iv (3)
      Integer :: ist, jst, ispn
      Integer :: i, j, k, l, m, n
      Integer :: i1, i2, i3, ir
      Real (8) :: vpc (3), vpql (3), vpqc (3)
      Real (8) :: v1 (3), v2 (3), v3 (3), t1
      Complex (8) zt1
! allocatable arrays
      Integer, Allocatable :: igpig (:)
      Integer, Allocatable :: igpqig (:)
      Real (8), Allocatable :: vgpl (:, :)
      Real (8), Allocatable :: vgpc (:, :)
      Real (8), Allocatable :: gpc (:)
      Real (8), Allocatable :: tpgpc (:, :)
      Real (8), Allocatable :: vgpql (:, :)
      Real (8), Allocatable :: vgpqc (:, :)
      Real (8), Allocatable :: gpqc (:)
      Real (8), Allocatable :: tpgpqc (:, :)
      Complex (8), Allocatable :: sfacgp (:, :)
      Complex (8), Allocatable :: sfacgpq (:, :)
      Complex (8), Allocatable :: apwalm1 (:, :, :, :)
      Complex (8), Allocatable :: apwalm2 (:, :, :, :)
      Complex (8), Allocatable :: evecfv1 (:, :)
      Complex (8), Allocatable :: evecfv2 (:, :)
      Complex (8), Allocatable :: evecsv1 (:, :)
      Complex (8), Allocatable :: evecsv2 (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Complex (8), Allocatable :: wfmt3 (:, :)
      Complex (8), Allocatable :: zfir1 (:)
      Complex (8), Allocatable :: zfir2 (:)
      Complex (8), Allocatable :: zfir3 (:)
      Complex (8), Allocatable :: zv (:)
      Complex (8), Allocatable :: epm (:, :, :)
! external functions
      Complex (8) zfmtinp, zdotc
      External zfmtinp, zdotc
      n = 3 * natmtot
! allocate local arrays
      Allocate (igpig(ngkmax))
      Allocate (igpqig(ngkmax))
      Allocate (vgpl(3, ngkmax))
      Allocate (vgpc(3, ngkmax))
      Allocate (gpc(ngkmax))
      Allocate (tpgpc(2, ngkmax))
      Allocate (vgpql(3, ngkmax))
      Allocate (vgpqc(3, ngkmax))
      Allocate (gpqc(ngkmax))
      Allocate (tpgpqc(2, ngkmax))
      Allocate (sfacgp(ngkmax, natmtot))
      Allocate (sfacgpq(ngkmax, natmtot))
      Allocate (apwalm1(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (apwalm2(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv1(nmatmax, nstfv))
      Allocate (evecfv2(nmatmax, nstfv))
      If (input%groundstate%tevecsv) Then
         Allocate (evecsv1(nstsv, nstsv))
         Allocate (evecsv2(nstsv, nstsv))
      End If
      Allocate (wfmt1(lmmaxapw, nrcmtmax))
      Allocate (wfmt2(lmmaxapw, nrcmtmax, nstfv))
      Allocate (wfmt3(lmmaxapw, nrcmtmax))
      Allocate (zfir1(ngrtot))
      Allocate (zfir2(ngrtot))
      Allocate (zfir3(ngrtot))
      Allocate (zv(ngkmax))
      Allocate (epm(nstfv, nstfv, n))
! p-vector in Cartesian coordinates
      Call r3mv (bvec, vpl, vpc)
! generate the G+p vectors
      Call gengpvec (vpl, vpc, ngp, igpig, vgpl, vgpc, gpc, tpgpc)
! generate the structure factors
      Call gensfacgp (ngp, vgpc, ngkmax, sfacgp)
! find the matching coefficients for k-point p
      Call match (ngp, gpc, tpgpc, sfacgp, apwalm1)
! get the eigenvectors for k-point p
      Call getevecfv (vpl, vgpl, evecfv1)
! p+q-vector in lattice coordinates
      vpql (:) = vpl (:) + vql (:, iq)
! map vector components to [0,1) interval
      Call r3frac (input%structure%epslat, vpql, iv)
! p+q-vector in Cartesian coordinates
      Call r3mv (bvec, vpql, vpqc)
! generate the G+p+q-vectors
      Call gengpvec (vpql, vpqc, ngpq, igpqig, vgpql, vgpqc, gpqc, &
     & tpgpqc)
! generate the structure factors
      Call gensfacgp (ngpq, vgpqc, ngkmax, sfacgpq)
! find the matching coefficients for k-point p+q
      Call match (ngpq, gpqc, tpgpqc, sfacgpq, apwalm2)
! get the eigenvectors for k-point p+q
      Call getevecfv (vpql, vgpql, evecfv2)
! set the first-variational matrix element array to zero
      epm (:, :, :) = 0.d0
!------------------------------------!
!     muffin-tin matrix elements     !
!------------------------------------!
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, nstfv
! calculate the wavefunction for k-point p+q
               Call wavefmt (input%groundstate%lradstep, &
              & input%groundstate%lmaxapw, is, ia, ngpq, apwalm2, &
              & evecfv2(:, ist), lmmaxapw, wfmt1)
! convert from spherical harmonics to spherical coordinates
               Call zgemm ('N', 'N', lmmaxapw, nrc, lmmaxapw, zone, &
              & zbshtapw, lmmaxapw, wfmt1, lmmaxapw, zzero, wfmt2(:, :, &
              & ist), lmmaxapw)
            End Do
            Do jst = 1, nstfv
! calculate the wavefunction for k-point p
               Call wavefmt (input%groundstate%lradstep, &
              & input%groundstate%lmaxapw, is, ia, ngp, apwalm1, &
              & evecfv1(:, jst), lmmaxapw, wfmt1)
! convert from spherical harmonics to spherical coordinates
               Call zgemm ('N', 'N', lmmaxapw, nrc, lmmaxapw, zone, &
              & zbshtapw, lmmaxapw, wfmt1, lmmaxapw, zzero, wfmt3, &
              & lmmaxapw)
! loop over phonon branches
               Do i = 1, n
! multiply the wavefunction by the change in effective potential
                  Do irc = 1, nrc
                     wfmt1 (:, irc) = wfmt3 (:, irc) * dveffmt (:, irc, &
                    & ias, i)
                  End Do
! add to the first-variational matrix elements
                  Do ist = 1, nstfv
                     epm (ist, jst, i) = epm (ist, jst, i) + zfmtinp &
                    & (.False., input%groundstate%lmaxapw, nrc, rcmt(:, &
                    & is), lmmaxapw, wfmt2(:, :, ist), wfmt1)
                  End Do
               End Do
            End Do
! end loops over atoms and species
         End Do
      End Do
!--------------------------------------!
!     interstitial matrix elements     !
!--------------------------------------!
! store G=q+p-p', where p' is the p+q-vector mapped to [0,1)
      v1 (:) = vqc (:, iq) + vpc (:) - vpqc (:)
! compute exp(i(q+p-p').r) for each r-vector on the grid
      ir = 0
      Do i3 = 0, ngrid (3) - 1
         v2 (3) = dble (i3) / dble (ngrid(3))
         Do i2 = 0, ngrid (2) - 1
            v2 (2) = dble (i2) / dble (ngrid(2))
            Do i1 = 0, ngrid (1) - 1
               v2 (1) = dble (i1) / dble (ngrid(1))
               ir = ir + 1
               Call r3mv (input%structure%crystal%basevect, v2, v3)
               t1 = v1 (1) * v3 (1) + v1 (2) * v3 (2) + v1 (3) * v3 (3)
               zfir1 (ir) = cmplx (Cos(t1), Sin(t1), 8)
            End Do
         End Do
      End Do
! compute interstitial wavefunctions for k-point p
      Do jst = 1, nstfv
         zfir2 (:) = 0.d0
         Do igp = 1, ngp
            ifg = igfft (igpig(igp))
            zfir2 (ifg) = evecfv1 (igp, jst)
         End Do
! Fourier transform wavefunction to real-space
         Call zfftifc (3, ngrid, 1, zfir2)
! multiply with the phase factor
         zfir2 (:) = zfir2 (:) * zfir1 (:)
! loop over phonon branches
         Do i = 1, n
! multiply the wavefunction with the change in effective potential
            zfir3 (:) = zfir2 (:) * dveffir (:, i)
! Fourier transform to G-space
            Call zfftifc (3, ngrid,-1, zfir3)
! store as wavefunction with G+p+q index
            Do igp = 1, ngpq
               ifg = igfft (igpqig(igp))
               zv (igp) = zfir3 (ifg)
            End Do
! add to the first-variational matrix elements
            Do ist = 1, nstfv
               epm (ist, jst, i) = epm (ist, jst, i) + zdotc (ngpq, &
              & evecfv2(:, ist), 1, zv, 1)
            End Do
         End Do
      End Do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
      If (input%groundstate%tevecsv) Then
! get the second-variational eigenvectors
         Call getevecsv (vpl, evecsv1)
         Call getevecsv (vpql, evecsv2)
         epmat (:, :, :) = 0.d0
         Do i = 1, nstsv
            Do j = 1, nstsv
               k = 0
               Do ispn = 1, nspinor
                  Do ist = 1, nstfv
                     k = k + 1
                     l = (ispn-1) * nstfv
                     Do jst = 1, nstfv
                        l = l + 1
                        zt1 = conjg (evecsv2(k, i)) * evecsv1 (l, j)
                        Do m = 1, n
                           epmat (i, j, m) = epmat (i, j, m) + epm &
                          & (ist, jst, m) * zt1
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      Else
         epmat (:, :, :) = epm (:, :, :)
      End If
      Deallocate (igpig, igpqig, vgpl, vgpc, gpc, tpgpc, vgpql, vgpqc, &
     & gpqc, tpgpqc)
      Deallocate (sfacgp, sfacgpq, apwalm1, apwalm2, evecfv1, evecfv2)
      If (input%groundstate%tevecsv) deallocate (evecsv1, evecsv2)
      Deallocate (wfmt1, wfmt2, wfmt3, zfir1, zfir2, zfir3, zv, epm)
      Return
End Subroutine
