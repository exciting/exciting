!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine epcouple
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ia, ias, ip
      Integer :: js, ja, jas
      Integer :: i, j, n, iv (3), isym
      Integer :: iq, ik, jk, ikq
      Integer :: ist, jst, ir, irc
      Real (8) :: vkql (3), x
      Real (8) :: t1, t2, t3, t4
      Complex (8) zt1
! allocatable arrays
      Real (8), Allocatable :: wq (:, :)
      Real (8), Allocatable :: gq (:, :)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: ev (:, :)
      Complex (8), Allocatable :: dveffmt (:, :, :, :)
      Complex (8), Allocatable :: dveffir (:, :)
      Complex (8), Allocatable :: zfmt (:, :, :)
      Complex (8), Allocatable :: gzfmt (:, :, :, :)
      Complex (8), Allocatable :: zflm (:)
      Complex (8), Allocatable :: zfir (:)
      Complex (8), Allocatable :: epmat (:, :, :)
      Complex (8), Allocatable :: gmq (:, :, :)
      Complex (8), Allocatable :: b (:, :)
! external functions
      Real (8) :: sdelta, gaunt
      External sdelta, gaunt
! initialise universal variables
      Call init0
      Call init1
      Call init2
! check k-point grid is commensurate with q-point grid
      iv (:) = Mod (input%groundstate%ngridk(:), ngridq(:))
      If ((iv(1) .Ne. 0) .Or. (iv(2) .Ne. 0) .Or. (iv(3) .Ne. 0)) Then
         Write (*,*)
         Write (*, '("Error(epcouple): k-point grid incommensurate with&
        & q-point grid")')
         Write (*, '(" ngridk : ", 3I6)') input%groundstate%ngridk
         Write (*, '(" ngridq : ", 3I6)') ngridq
         Write (*,*)
         Stop
      End If
      n = 3 * natmtot
! allocate local arrays
      Allocate (wq(n, nqpt))
      Allocate (gq(n, nqpt))
      Allocate (dynq(n, n, nqpt))
      Allocate (ev(n, n))
      Allocate (dveffmt(lmmaxapw, nrcmtmax, natmtot, n))
      Allocate (dveffir(ngrtot, n))
      Allocate (zfmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (gzfmt(lmmaxvr, nrcmtmax, 3, natmtot))
      Allocate (zflm(lmmaxapw))
      Allocate (zfir(ngrtot))
      Allocate (gmq(n, n, nqpt))
      Allocate (b(n, n))
! read in the density and potentials from file
      Call readstate
! read Fermi energy from file
      Call readfermi
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! get the eigenvalues from file
      Do ik = 1, nkpt
         Call getevalsv (vkl(:, ik), evalsv(:, ik))
      End Do
! compute the occupancies and density of states at the Fermi energy
      Call occupy
! read in the dynamical matrices
      Call readdyn (.true.,dynq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! compute the gradients of the effective potential for the rigid-ion term
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! convert potential to complex spherical harmonic expansion
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call rtozflm (input%groundstate%lmaxvr, veffmt(:, ir, &
              & ias), zfmt(:, irc, ias))
            End Do
            Call gradzfmt (input%groundstate%lmaxvr, nrcmt(is), rcmt(:, &
           & is), lmmaxvr, nrcmtmax, zfmt(:, :, ias), gzfmt(:, :, :, &
           & ias))
         End Do
      End Do
! loop over phonon q-points
      Do iq = 1, nqpt
         Write (*, '("Info(epcouple): ", I6, " of ", I6, " q-points")') &
        & iq, nqpt
! diagonalise the dynamical matrix
         Call dyndiag (dynq(:, :, iq), wq(:, iq), ev)
! loop over phonon branches
         Do j = 1, n
! find change effective potential for mode j
            dveffmt (:, :, :, j) = 0.d0
            dveffir (:, j) = 0.d0
            i = 0
            Do is = 1, nspecies
! prefactor
               t1 = 2.d0 * spmass (is) * wq (j, iq)
               If (t1 .Gt. 1.d-8) Then
                  t1 = 1.d0 / Sqrt (t1)
               Else
                  t1 = 0.d0
               End If
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  Do ip = 1, 3
                     i = i + 1
! read in the Cartesian change in effective potential
                     Call readdveff (iq, is, ia, ip, zfmt, zfir)
! add the rigid-ion term
                     Do irc = 1, nrcmt (is)
                        zfmt (:, irc, ias) = zfmt (:, irc, ias) - gzfmt &
                       & (:, irc, ip, ias)
                     End Do
! multiply with eigenvector component and add to total
                     zt1 = t1 * ev (i, j)
                     Do js = 1, nspecies
                        Do ja = 1, natoms (js)
                           jas = idxas (ja, js)
                           Do irc = 1, nrcmt (js)
                              dveffmt (1:lmmaxvr, irc, jas, j) = &
                             & dveffmt (1:lmmaxvr, irc, jas, j) + zt1 * &
                             & zfmt (1:lmmaxvr, irc, jas)
                           End Do
                        End Do
                     End Do
                     dveffir (:, j) = dveffir (:, j) + zt1 * zfir (:)
                  End Do
               End Do
            End Do
! multiply the interstitial potential with the characteristic function
            dveffir (:, j) = dveffir (:, j) * cfunir (:)
! convert muffin-tin potential to spherical coordinates on the lmaxapw covering
            zflm (:) = 0.d0
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  Do irc = 1, nrcmt (is)
                     zflm (1:lmmaxvr) = dveffmt (1:lmmaxvr, irc, ias, &
                    & j)
                     Call zgemv ('N', lmmaxapw, lmmaxapw, zone, &
                    & zbshtapw, lmmaxapw, zflm, 1, zzero, dveffmt(:, &
                    & irc, ias, j), 1)
                  End Do
               End Do
            End Do
         End Do
! zero the phonon linewidths array
         gq (:, iq) = 0.d0
! begin parallel loop over non-reduced k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(epmat,jk,vkql,isym) &
!$OMP PRIVATE(ikq,ist,jst,i) &
!$OMP PRIVATE(x,t1,t2,t3,t4)
!$OMP DO
         Do ik = 1, nkptnr
            Allocate (epmat(nstsv, nstsv, n))
! equivalent reduced k-point
            jk = ikmap (ivknr(1, ik), ivknr(2, ik), ivknr(3, ik))
! compute the electron-phonon coupling matrix elements
            Call genepmat (iq, vklnr(:, ik), dveffmt, dveffir, epmat)
! k+q-vector in lattice coordinates
            vkql (:) = vklnr (:, ik) + vql (:, iq)
! index to k+q-vector
            Call findkpt (vkql, isym, ikq)
            t1 = twopi * wkptnr (ik) * (occmax/2.d0)
! loop over second-variational states
            Do ist = 1, nstsv
               x = (evalsv(ist, ikq)-efermi) / input%groundstate%swidth
               t2 = sdelta (input%groundstate%stypenumber, x) / &
              & input%groundstate%swidth
! loop over phonon branches
               Do i = 1, n
                  Do jst = 1, nstsv
                     x = (evalsv(jst, jk)-efermi) / &
                    & input%groundstate%swidth
                     t3 = sdelta (input%groundstate%stypenumber, x) / &
                    & input%groundstate%swidth
                     t4 = dble (epmat(ist, jst, i)) ** 2 + aimag &
                    & (epmat(ist, jst, i)) ** 2
!$OMP CRITICAL
                     gq (i, iq) = gq (i, iq) + wq (i, iq) * t1 * t2 * &
                    & t3 * t4
!$OMP END CRITICAL
                  End Do
               End Do
            End Do
            Deallocate (epmat)
! end loop over k-points
         End Do
!$OMP END DO
!$OMP END PARALLEL
! end loop over phonon q-points
      End Do
! write the phonon linewidths to file
      Call writegamma (gq)
! write electron-phonon coupling constants to file
      Call writelambda (wq, gq)
      Deallocate (wq, gq, dynq, ev, dveffmt, dveffir)
      Deallocate (zfmt, gzfmt, zflm, zfir, gmq, b)
      Return
End Subroutine
