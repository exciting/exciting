!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genpmat
! !INTERFACE:
!
!
Subroutine genpmat (ngp, igpig, vgpc, apwalm, evecfv, evecsv, pmat)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmat   : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!   Fixed bug found by Juergen Spitaler, September 2006 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv)
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
      Complex (8), Intent (Out) :: pmat (3, nstsv, nstsv)
! local variables
      Integer :: ispn, is, ia, ist, jst
      Integer :: i, j, k, l, igp, ifg, ir
      Complex (8) zsum, zt1, zv (3)
! allocatable arrays
      Complex (8), Allocatable :: wfmt (:, :, :)
      Complex (8), Allocatable :: gwfmt (:, :, :, :)
      Complex (8), Allocatable :: wfir (:, :)
      Complex (8), Allocatable :: gwfir (:, :, :)
      Complex (8), Allocatable :: pm (:, :, :)
! external functions
      Complex (8) zfmtinp
      External zfmtinp
      Allocate (wfmt(lmmaxapw, nrcmtmax, nstfv))
      Allocate (gwfmt(lmmaxapw, nrcmtmax, 3, nstfv))
      Allocate (wfir(ngrtot, nstfv))
      Allocate (gwfir(ngrtot, 3, nstfv))
      Allocate (pm(3, nstfv, nstfv))
! set the momentum matrix elements to zero
      pm (:, :, :) = 0.d0
! calculate momentum matrix elements in the muffin-tin
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            Do ist = 1, nstfv
! calculate the wavefunction
               Call wavefmt (input%groundstate%lradstep, &
              & input%groundstate%lmaxapw, is, ia, ngp, apwalm, &
              & evecfv(:, ist), lmmaxapw, wfmt(:, :, ist))
! calculate the gradient
               Call gradzfmt (input%groundstate%lmaxapw, nrcmt(is), &
              & rcmt(:, is), lmmaxapw, nrcmtmax, wfmt(:, :, ist), &
              & gwfmt(:, :, :, ist))
            End Do
            Do ist = 1, nstfv
               Do jst = ist, nstfv
                  Do i = 1, 3
                     zt1 = zfmtinp (.True., input%groundstate%lmaxapw, &
                    & nrcmt(is), rcmt(:, is), lmmaxapw, wfmt(:, :, &
                    & ist), gwfmt(:, :, i, jst))
                     pm (i, ist, jst) = pm (i, ist, jst) + zt1
                  End Do
               End Do
            End Do
         End Do
      End Do
! calculate momemntum matrix elements in the interstitial region
      wfir (:, :) = 0.d0
      gwfir (:, :, :) = 0.d0
      Do ist = 1, nstfv
         Do igp = 1, ngp
            ifg = igfft (igpig(igp))
            zt1 = evecfv (igp, ist)
            wfir (ifg, ist) = zt1
! calculate the gradient
            Do i = 1, 3
               gwfir (ifg, i, ist) = zi * vgpc (i, igp) * zt1
            End Do
         End Do
! Fourier transform the wavefunction to real-space
         Call zfftifc (3, ngrid, 1, wfir(:, ist))
         Do i = 1, 3
            Call zfftifc (3, ngrid, 1, gwfir(:, i, ist))
         End Do
      End Do
! find the overlaps
      Do ist = 1, nstfv
         Do jst = ist, nstfv
            Do i = 1, 3
               zsum = 0.d0
               Do ir = 1, ngrtot
                  zsum = zsum + cfunir (ir) * conjg (wfir(ir, ist)) * &
                 & gwfir (ir, i, jst)
               End Do
               zt1 = zsum / dble (ngrtot)
               pm (i, ist, jst) = pm (i, ist, jst) + zt1
            End Do
         End Do
      End Do
! multiply by -i and set lower triangular part
      Do ist = 1, nstfv
         Do jst = ist, nstfv
            pm (:, ist, jst) = - zi * pm (:, ist, jst)
            pm (:, jst, ist) = conjg (pm(:, ist, jst))
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
                        zv (:) = zv (:) + zt1 * pm (:, ist, jst)
                     End Do
                  End Do
               End Do
               pmat (:, i, j) = zv (:)
            End Do
         End Do
      Else
         pmat (:, :, :) = pm (:, :, :)
      End If
      Deallocate (wfmt, gwfmt, wfir, gwfir, pm)
      Return
End Subroutine
!EOC
