!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genshtmat
! !INTERFACE:
!
!
Subroutine genshtmat
! !USES:
      Use modinput
      Use mod_muffin_tin, only: lmmaxapw, lmmaxvr
      Use mod_SHT, only: rfshtapw, rbshtapw, rfshtvr, rbshtvr, &
                         zfshtapw, zbshtapw, zfshtvr, zbshtvr, &
                         gen_rshtmat, gen_zshtmat
#ifdef XS
      Use modxs
#endif
! !DESCRIPTION:
!   Generates the forward and backward spherical harmonic transformation (SHT)
!   matrices using the spherical covering set produced by the routine
!   {\tt sphcover}. These matrices are used to transform a function between its
!   $(l,m)$-expansion coefficients and its values at the $(\theta,\phi)$ points
!   on the sphere.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Modified July 2022 (SeTi)
!EOP
!BOC
      Implicit None
! local variables
! allocatable arrays
      Real (8), Allocatable :: tp (:, :)

! generate matrices for lmaxapw
      call gen_rshtmat( input%groundstate%lmaxapw, rfshtapw, rbshtapw)
      call gen_zshtmat( input%groundstate%lmaxapw, zfshtapw, zbshtapw)
! generate matrices for lmaxvr
      call gen_rshtmat( input%groundstate%lmaxvr, rfshtvr, rbshtvr)
      call gen_zshtmat( input%groundstate%lmaxvr, zfshtvr, zbshtvr)
! generate spherical covering set for lmaxvr
      Allocate (tp(2, lmmaxapw))
      Call sphcover (lmmaxvr, tp)
#ifdef XS
      If (allocated(sphcov)) deallocate (sphcov)
      Allocate (sphcov(3, lmmaxapw))
      If (allocated(sphcovtp)) deallocate (sphcovtp)
      Allocate (sphcovtp(2, lmmaxapw))
      sphcovtp (:, :) = tp (:, :)
      sphcov (1, :) = Sin (tp(1, :)) * Cos (tp(2, :))
      sphcov (2, :) = Sin (tp(1, :)) * Sin (tp(2, :))
      sphcov (3, :) = Cos (tp(1, :))
#endif
      Deallocate (tp)

End Subroutine
!EOC
