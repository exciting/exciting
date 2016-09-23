
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: chgdist
! !INTERFACE:
Subroutine chgdist(rhomtref,rhoirref)
! !USES:
use modmain
! !DESCRIPTION:
!   Calculated the charge distance between two charge densities of the current
!   and of the last iteration according to
!   the expression
!   $$
!     \Delta Q = \int_{\Omega} {\rm d}^3 r [\rho^{(n)}({\bf r}) -
!     \rho^{(n-1)}({\bf r})].
!   $$
!   Based on the routine {\tt charge}.
!
! !REVISION HISTORY:
!   Created 2010 (Sagmeister)
!   Modified 2014 (DIN)
!EOP
!BOC
      implicit none
      real(8), intent(in) :: rhomtref(lmmaxvr,nrmtmax,natmtot)
      real(8), intent(in) :: rhoirref(ngrtot)
! local variables
      Integer :: is, ia, ias, ir
      real(8) :: sum, chgdstmt, chgdstir
! automatic arrays
      real(8) :: fr (nrmtmax), gr (nrmtmax), hr1(lmmaxvr), hr2(lmmaxvr)
      real(8) :: cf (3, nrmtmax), sht00(lmmaxvr)
! external functions
      real(8), external :: ddot
! find the muffin-tin charges
      chgdstmt=0.d0
      chgdstir=0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               hr1(:)=rhomtref(:,ir,ias)-rhomt(:,ir,ias)
               ! backward transform to real space
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0,rbshtvr, lmmaxvr, &
                 hr1, 1,0.d0, hr2, 1)
               ! take absolute value
               hr2(:)=abs(hr2(:))
               ! forward transform for 00 component
               sht00(:)=rfshtvr(1,:)
               fr(ir)=ddot(lmmaxvr,hr2,1,sht00,1) * spr (ir, is) ** 2
            End Do
            Call fderiv (-1, nrmt(is), spr(:, is), fr, gr, cf)
            chgdstmt = chgdstmt + fourpi * y00 * gr (nrmt(is))
         End Do
      End Do
! find the interstitial charge
      sum = 0.d0
      Do ir = 1, ngrtot
         sum = sum + abs(rhoirref(ir) - rhoir(ir)) * cfunir(ir)
      End Do
      chgdstir = sum * omega / dble (ngrtot)
      chgdst = (chgdstmt + chgdstir)/chgtot
end subroutine
!EOC
