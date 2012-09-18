!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
!
!
Subroutine rhoinit
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
! polynomial order of smooth step function
      Integer, Parameter :: n = 4
      Integer :: lmax, lmmax, l, m, lm, ir, irc
      Integer :: is, ia, ias, ig, ifg
      Real (8) :: x, t1, t2
      Complex (8) zt1, zt2, zt3
! automatic arrays
      Real (8) :: fr (spnrmax), gr (spnrmax), cf (3, spnrmax)
! allocatable arrays
      Real (8), Allocatable :: jlgr (:, :)
      Real (8), Allocatable :: th (:, :)
      Real (8), Allocatable :: ffacg (:)
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: zfft (:)
! maximum angular momentum for density initialisation
      lmax = 1
      lmmax = (lmax+1) ** 2
! allocate local arrays
      Allocate (jlgr(0:lmax, nrcmtmax))
      Allocate (th(nrmtmax, nspecies))
      Allocate (ffacg(ngvec))
      Allocate (zfmt(lmmax, nrcmtmax))
      Allocate (zfft(ngrtot))
! zero the charge density and magnetisation arrays
      rhomt (:, :, :) = 0.d0
      rhoir (:) = 0.d0
      If (associated(input%groundstate%spin)) Then
         magmt (:, :, :, :) = 0.d0
         magir (:, :) = 0.d0
      End If
! compute the superposition of all the atomic density tails
      zfft (:) = 0.d0
      Do is = 1, nspecies
! generate smooth step function
         Do ir = 1, nrmt (is)
            th (ir, is) = (spr(ir, is)/rmt(is)) ** n
         End Do
         Do ig = 1, ngvec
            Do ir = 1, spnr (is)
               x = gc (ig) * spr (ir, is)
               Call sbessel (0, x, jlgr(0, 1))
               t1 = sprho (ir, is) * jlgr (0, 1) * spr (ir, is) ** 2
               If (ir .Lt. nrmt(is)) Then
                  fr (ir) = th (ir, is) * t1
               Else
                  fr (ir) = t1
               End If
            End Do
            Call fderiv (-1, spnr(is), spr(:, is), fr, gr, cf)
            ffacg (ig) = (fourpi/omega) * gr (spnr(is))
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ig = 1, ngvec
               ifg = igfft (ig)
               zfft (ifg) = zfft (ifg) + ffacg (ig) * conjg (sfacg(ig, &
              & ias))
            End Do
         End Do
      End Do
! compute the tails in each muffin-tin
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zfmt (:, :) = 0.d0
            Do ig = 1, ngvec
               ifg = igfft (ig)
               Do irc = 1, nrcmt (is)
                  x = gc (ig) * rcmt (irc, is)
                  Call sbessel (lmax, x, jlgr(:, irc))
               End Do
               zt1 = fourpi * zfft (ifg) * sfacg (ig, ias)
               lm = 0
               Do l = 0, lmax
                  zt2 = zt1 * zil (l)
                  Do m = - l, l
                     lm = lm + 1
                     zt3 = zt2 * conjg (ylmg(lm, ig))
                     Do irc = 1, nrcmt (is)
                        zfmt (lm, irc) = zfmt (lm, irc) + jlgr (l, irc) &
                       & * zt3
                     End Do
                  End Do
               End Do
            End Do
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call ztorflm (lmax, zfmt(:, irc), rhomt(:, ir, ias))
            End Do
         End Do
      End Do
! convert the density from a coarse to a fine radial mesh
      Call rfmtctof (rhomt)
! add the atomic charge density and the excess charge in each muffin-tin
      t1 = input%groundstate%chgexs / omega
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               t2 = (t1+(1.d0-th(ir, is))*sprho(ir, is)) / y00
               rhomt (1, ir, ias) = rhomt (1, ir, ias) + t2
            End Do
         End Do
      End Do
! interstitial density determined from the atomic tails and excess charge
      Call zfftifc (3, ngrid, 1, zfft)
      Do ir = 1, ngrtot
         rhoir (ir) = dble (zfft(ir)) + t1
      End Do
! compute the total charge
      Call charge
! normalise the density
      Call rhonorm
      Deallocate (jlgr, th, ffacg, zfmt, zfft)
      Return
End Subroutine
!EOC
