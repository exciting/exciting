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
      Use modmain
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
!EOP
!BOC
      Implicit None
! local variables
      Integer :: itp, lwork, info
! allocatable arrays
      Integer, Allocatable :: ipiv (:)
      Real (8), Allocatable :: tp (:, :)
      Real (8), Allocatable :: rlm (:)
      Real (8), Allocatable :: work (:)
      Complex (8), Allocatable :: ylm (:)
      Complex (8), Allocatable :: zwork (:)
      Allocate (tp(2, lmmaxapw))
      Allocate (rlm(lmmaxapw))
      Allocate (ylm(lmmaxapw))
      Allocate (ipiv(lmmaxapw))
      lwork = 2 * lmmaxapw
      Allocate (work(lwork))
      Allocate (zwork(lwork))
! allocate real SHT matrices for lmaxapw
      If (allocated(rbshtapw)) deallocate (rbshtapw)
      Allocate (rbshtapw(lmmaxapw, lmmaxapw))
      If (allocated(rfshtapw)) deallocate (rfshtapw)
      Allocate (rfshtapw(lmmaxapw, lmmaxapw))
! allocate real SHT matrices for lmaxvr
      If (allocated(rbshtvr)) deallocate (rbshtvr)
      Allocate (rbshtvr(lmmaxvr, lmmaxvr))
      If (allocated(rfshtvr)) deallocate (rfshtvr)
      Allocate (rfshtvr(lmmaxvr, lmmaxvr))
! allocate complex SHT matrices for lmaxapw
      If (allocated(zbshtapw)) deallocate (zbshtapw)
      Allocate (zbshtapw(lmmaxapw, lmmaxapw))
      If (allocated(zfshtapw)) deallocate (zfshtapw)
      Allocate (zfshtapw(lmmaxapw, lmmaxapw))
! allocate complex SHT matrices for lmaxvr
      If (allocated(zbshtvr)) deallocate (zbshtvr)
      Allocate (zbshtvr(lmmaxvr, lmmaxvr))
      If (allocated(zfshtvr)) deallocate (zfshtvr)
      Allocate (zfshtvr(lmmaxvr, lmmaxvr))
! generate spherical covering set for lmaxapw
      Call sphcover (lmmaxapw, tp)
! generate real and complex spherical harmonics and set the backward SHT arrays
      Do itp = 1, lmmaxapw
         Call genrlm (input%groundstate%lmaxapw, tp(:, itp), rlm)
         rbshtapw (itp, 1:lmmaxapw) = rlm (1:lmmaxapw)
         Call genylm (input%groundstate%lmaxapw, tp(:, itp), ylm)
         zbshtapw (itp, 1:lmmaxapw) = ylm (1:lmmaxapw)
      End Do
! generate spherical covering set for lmaxvr
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
! generate real and complex spherical harmonics and set the backward SHT arrays
      Do itp = 1, lmmaxvr
         Call genrlm (input%groundstate%lmaxvr, tp(:, itp), rlm)
         rbshtvr (itp, 1:lmmaxvr) = rlm (1:lmmaxvr)
         Call genylm (input%groundstate%lmaxvr, tp(:, itp), ylm)
         zbshtvr (itp, 1:lmmaxvr) = ylm (1:lmmaxvr)
      End Do
! find the forward SHT arrays
! real, lmaxapw
      rfshtapw (:, :) = rbshtapw (:, :)
      Call dgetrf (lmmaxapw, lmmaxapw, rfshtapw, lmmaxapw, ipiv, info)
      If (info .Ne. 0) Go To 10
      Call dgetri (lmmaxapw, rfshtapw, lmmaxapw, ipiv, work, lwork, &
     & info)
      If (info .Ne. 0) Go To 10
! real, lmaxvr
      rfshtvr (:, :) = rbshtvr (:, :)
      Call dgetrf (lmmaxvr, lmmaxvr, rfshtvr, lmmaxvr, ipiv, info)
      If (info .Ne. 0) Go To 10
      Call dgetri (lmmaxvr, rfshtvr, lmmaxvr, ipiv, work, lwork, info)
      If (info .Ne. 0) Go To 10
! complex, lmaxapw
      zfshtapw (:, :) = zbshtapw (:, :)
      Call zgetrf (lmmaxapw, lmmaxapw, zfshtapw, lmmaxapw, ipiv, info)
      If (info .Ne. 0) Go To 10
      Call zgetri (lmmaxapw, zfshtapw, lmmaxapw, ipiv, zwork, lwork, &
     & info)
      If (info .Ne. 0) Go To 10
! complex, lmaxvr
      zfshtvr (:, :) = zbshtvr (:, :)
      Call zgetrf (lmmaxvr, lmmaxvr, zfshtvr, lmmaxvr, ipiv, info)
      If (info .Ne. 0) Go To 10
      Call zgetri (lmmaxvr, zfshtvr, lmmaxvr, ipiv, zwork, lwork, info)
      If (info .Ne. 0) Go To 10
      Deallocate (tp, rlm, ylm, ipiv, work, zwork)
      Return
10    Continue
      Write (*,*)
      Write (*, '("Error(genshtmat): unable to find inverse spherical h&
     &armonic transform")')
      Write (*, '(" => improper spherical covering")')
      Write (*,*)
      Stop
End Subroutine
!EOC
