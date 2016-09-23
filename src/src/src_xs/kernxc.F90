!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: kernxc
! !INTERFACE:
!
!
Subroutine kernxc
! !DESCRIPTION:
!   Computes the ALDA exchange-correlation kernel. In the
!   muffin-tin, the density is transformed from spherical harmonic coefficients
!   $\rho_{lm}$ to spherical coordinates $(\theta,\phi)$ with a backward
!   spherical harmonic transformation (SHT). Once calculated, the
!   exchange-correlation potential and energy density are transformed with a
!   forward SHT. This routine is based upon the routine {\tt potxc}.
!
! !REVISION HISTORY:
!   Created March 2007 (Sagmeister)
!EOP
!BOC
      Use modmain
      Use modxs
      Implicit None
  ! local variables
      Real (8), Allocatable :: dvx (:, :), dvc (:, :), rftp (:, :)
      Real (8), Allocatable :: f1ir (:), f2ir (:), f1mt (:, :, :), f2mt &
     & (:, :, :)
      Complex (8), Allocatable :: zftp (:, :)
      Integer :: m, is, ia, ias, ir, itp
      Allocate (rftp(lmmaxvr, 1), zftp(lmmaxvr, 1))
      m = Max (lmmaxvr, ngrtot)
      Allocate (dvx(m, 1), dvc(m, 1))
  !-----------------------------!
  !     muffin-tin potential    !
  !-----------------------------!
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
           !--------------------------!
           !     spin-unpolarised     !
           !--------------------------!
           ! convert density to real space
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, rhomt(1, ir, ias), 1, 0.d0, rftp, 1)
               Call xcd_pwca (lmmaxvr, rftp, dvx, dvc)
               Do itp = 1, lmmaxvr
                  rftp (itp, 1) = dvx (itp, 1) + dvc (itp, 1)
               End Do
           ! convert kernel to spherical-harmonics expansion
               zftp (:, :) = rftp (:, :)
               Call zgemv ('N', lmmaxvr, lmmaxvr, zone, zfshtvr, &
              & lmmaxvr, zftp, 1, zzero, fxcmt(1, ir, ias), 1)
            End Do
         End Do
      End Do
  !--------------------------------!
  !     interstitial potential     !
  !--------------------------------!
  !--------------------------!
  !     spin-unpolarised     !
  !--------------------------!
      Call xcd_pwca (ngrtot, rhoir, dvx, dvc)
      fxcir (1:ngrtot) = dvx (1:ngrtot, 1) + dvc (1:ngrtot, 1)
      Allocate (f1ir(ngrtot), f1mt(lmmaxvr, nrmtmax, natmtot))
      Allocate (f2ir(ngrtot), f2mt(lmmaxvr, nrmtmax, natmtot))
      f1ir (:) = dble (fxcir(:))
      f1mt (:, :, :) = dble (fxcmt(:, :, :))
      f2ir (:) = aimag (fxcir(:))
      f2mt (:, :, :) = aimag (fxcmt(:, :, :))
  ! symmetrise the exchange-correlation kernel
      Call symrf (1, f1mt, f1ir)
      Call symrf (1, f2mt, f2ir)
  ! back-substitute
      fxcir (:) = f1ir (:) + zi * f2ir (:)
      fxcmt (:, :, :) = f1mt (:, :, :) + zi * f2mt (:, :, :)
      Deallocate (f1ir, f2ir, f1mt, f2mt)
      Deallocate (rftp, zftp, dvx, dvc)
End Subroutine kernxc
!EOC
