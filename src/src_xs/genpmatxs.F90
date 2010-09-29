!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genpmatxs
! !INTERFACE:
!
!
Subroutine genpmatxs (ngp, igpig, vgpc, evecfv, evecsv, pmat)
! !USES:
      Use modinput
      Use modmain
      Use modxs, Only: apwcmt, locmt, ripaa, ripalo, riploa, riplolo
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmat   : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!   The gradient is applied explicitly only to the radial functions and
!   corresponding spherical harmonics for the muffin-tin part. In the
!   interstitial region the gradient is evaluated analytically.
!   Parts taken from the routine {\tt genpmat}.
!
! !REVISION HISTORY:
!   Created April 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv)
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
      Complex (8), Intent (Out) :: pmat (3, nstsv, nstsv)
  ! local variables
      Integer :: ispn, is, ia, ias, ist, jst
      Integer :: ist1, l1, m1, lm1, l3, m3, lm3, io, io1, io2, ilo, &
     & ilo1, ilo2
      Integer :: i, j, k, l
      Integer :: igp1, igp2, ig1, ig2, ig, iv1 (3), iv (3)
      Complex (8) :: zt1, zv (3)
  ! allocatable arrays
      Complex (8), Allocatable :: wfmt (:, :, :)
      Complex (8), Allocatable :: gwfmt (:, :, :, :)
      Complex (8), Allocatable :: pm (:, :, :)
      Complex (8), Allocatable :: cfunt (:, :), h (:, :), pmt (:, :)
      Complex (8), Allocatable :: evecfv1 (:, :), evecfv2 (:, :)
      Complex (8), Allocatable :: zv2 (:)
  ! external functions
      Complex (8) zfmtinp
      External zfmtinp
      Allocate (zv2(nstfv))
      Allocate (wfmt(lmmaxapw, nrcmtmax, nstfv))
      Allocate (gwfmt(lmmaxapw, nrcmtmax, 3, nstfv))
      Allocate (cfunt(ngp, ngp))
      Allocate (h(ngp, nstfv))
      Allocate (pmt(nstfv, nstfv))
      Allocate (evecfv1(nstfv, ngp), evecfv2(ngp, nstfv))
      Allocate (pm(nstfv, nstfv, 3))
  ! set the momentum matrix elements to zero
      pm (:, :, :) = 0.d0
  ! loop over species and atoms
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
            Do j = 1, 3
               Do l1 = 0, input%groundstate%lmaxapw
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     Do io1 = 1, apword (l1, is)
                        zv2 (:) = zzero
                        Do l3 = 0, input%groundstate%lmaxapw
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Do io2 = 1, apword (l3, is)
                                 Call zaxpy (nstfv, zone*ripaa(io1, &
                                & lm1, io2, lm3, ias, j), apwcmt(1, &
                                & io2, lm3, ias), 1, zv2, 1)
                              End Do
                           End Do
                        End Do
                        Call zoutpr (nstfv, nstfv, zone, apwcmt(1, io1, &
                       & lm1, ias), zv2, pm(1, 1, j))
                     End Do
                  End Do
               End Do
            End Do
            If (nlotot .Gt. 0) Then
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
               Do j = 1, 3
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        Do io = 1, apword (l3, is)
                           zv2 (:) = zzero
                           Do ilo = 1, nlorb (is)
                              l1 = lorbl (ilo, is)
                              Do m1 = - l1, l1
                                 lm1 = idxlm (l1, m1)
                                 Call zaxpy (nstfv, zone*ripalo(io, &
                                & lm3, ilo, m1, ias, j), locmt(1, ilo, &
                                & m1, ias), 1, zv2, 1)
                              End Do
                           End Do
                           Call zoutpr (nstfv, nstfv, zone, apwcmt(1, &
                          & io, lm3, ias), zv2, pm(1, 1, j))
                        End Do
                     End Do
                  End Do
               End Do
           !--------------------------------------!
           !     local-orbital-APW contribution   !
           !--------------------------------------!
               Do j = 1, 3
                  Do ilo = 1, nlorb (is)
                     l1 = lorbl (ilo, is)
                     Do m1 = - l1, l1
                        lm1 = idxlm (l1, m1)
                        zv2 (:) = zzero
                        Do l3 = 0, input%groundstate%lmaxapw
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Do io = 1, apword (l3, is)
                                 Call zaxpy (nstfv, zone*riploa(ilo, &
                                & m1, io, lm3, ias, j), apwcmt(1, io, &
                                & lm3, ias), 1, zv2, 1)
                              End Do
                           End Do
                        End Do
                        Call zoutpr (nstfv, nstfv, zone, locmt(1, ilo, &
                       & m1, ias), zv2, pm(1, 1, j))
                     End Do
                  End Do
               End Do
           !------------------------------------------------!
           !     local-orbital-local-orbital contribution   !
           !------------------------------------------------!
               Do j = 1, 3
                  Do ilo1 = 1, nlorb (is)
                     l1 = lorbl (ilo1, is)
                     Do m1 = - l1, l1
                        lm1 = idxlm (l1, m1)
                        zv2 (:) = zzero
                        Do ilo2 = 1, nlorb (is)
                           l3 = lorbl (ilo2, is)
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Call zaxpy (nstfv, zone*riplolo(ilo1, m1, &
                             & ilo2, m3, ias, j), locmt(1, ilo2, m3, &
                             & ias), 1, zv2, 1)
                           End Do
                        End Do
                        Call zoutpr (nstfv, nstfv, zone, locmt(1, ilo1, &
                       & m1, ias), zv2, pm(1, 1, j))
                     End Do
                  End Do
               End Do
           ! end case of local orbitals
            End If
        ! end loop over atoms and species
         End Do
      End Do
  ! multiply y-component with imaginary unit
      pm (:, :, 2) = zi * pm (:, :, 2)
  !  calculate momentum matrix elements in the interstitial region
      Forall (ist1=1:nstfv)
         evecfv1 (ist1, :) = conjg (evecfv(1:ngp, ist1))
      End Forall
      evecfv2 (:, :) = evecfv (1:ngp, :)
      Do j = 1, 3
         Do igp1 = 1, ngp
            ig1 = igpig (igp1)
            iv1 (:) = ivg (:, ig1)
            Do igp2 = 1, ngp
               ig2 = igpig (igp2)
               iv (:) = iv1 (:) - ivg (:, ig2)
               ig = ivgig (iv(1), iv(2), iv(3))
               cfunt (igp1, igp2) = zi * vgpc (j, igp2) * cfunig (ig)
            End Do
         End Do
         Call zgemm ('n', 'n', ngp, nstfv, ngp, zone, cfunt, ngp, &
        & evecfv2, ngp, zzero, h, ngp)
         Call zgemm ('n', 'n', nstfv, nstfv, ngp, zone, evecfv1, nstfv, &
        & h, ngp, zzero, pmt, nstfv)
         pm (:, :, j) = pm (:, :, j) + pmt (:, :)
      End Do
  ! multiply by -i and set lower triangular part
      Do ist = 1, nstfv
         Do jst = ist, nstfv
            pm (ist, jst, :) = - zi * pm (ist, jst, :)
            pm (jst, ist, :) = conjg (pm(ist, jst, :))
         End Do
      End Do
  ! compute the second-variational momentum matrix elements
      If (input%groundstate%tevecsv) Then
         Do i = 1, nstsv
            Do j = 1, nstsv
               zv (:) = 0.d0
               k = 0
               Do ispn = 1, nspinor
                  Do ist = 1, nstfv
                     k = k + 1
                     l = (ispn-1) * nstfv
                     Do jst = 1, nstfv
                        l = l + 1
                        zt1 = conjg (evecsv(k, i)) * evecsv (l, j)
                        zv (:) = zv (:) + zt1 * pm (ist, jst, :)
                     End Do
                  End Do
               End Do
               pmat (:, i, j) = zv (:)
            End Do
         End Do
      Else
         Do j = 1, 3
            pmat (j, :, :) = pm (:, :, j)
         End Do
      End If
      Deallocate (wfmt, gwfmt, pm, cfunt, h, pmt, evecfv1, evecfv2)
End Subroutine genpmatxs
!EOC
