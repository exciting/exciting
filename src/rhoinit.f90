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
! zero the charge density and magnetisation arrays
      rhomt (:, :, :) = 0.d0
      rhoir (:) = 0.d0
      If (associated(input%groundstate%spin)) Then
         magmt (:, :, :, :) = 0.d0
         magir (:, :) = 0.d0
      End If
      
      t1=0.d0
      t2=0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            do ir = 1, nrmt(is)
              rhomt(1,ir,ias)=sprho(ir,is)/y00
              fr(ir)=sprho(ir,is)*(spr(ir,is)**2)
            enddo
            Call fderiv (-1, nrmt(is), spr(:, is), fr, gr, cf)
            t1 = t1 + fourpi * gr (nrmt(is))
            t2 = t2 + (fourpi/3)*(rmt(is)**3)
         enddo
      enddo
      do ir=1,ngrtot
        rhoir(ir)=(chgtot-t1)/(omega-t2)
      enddo

! compute the total charge
      Call charge
! normalise the density
      Call rhonorm
      Return
End Subroutine
!EOC
