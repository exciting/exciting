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
Subroutine genshtmat3
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
      Integer :: itp, lwork, info, i
      Complex (8) :: ss
! allocatable arrays
      Integer, Allocatable :: ipiv (:)
      Real (8), Allocatable :: tp (:, :)
      Real (8), Allocatable :: rlm (:)
      Real (8), Allocatable :: work (:)
      Complex (8), Allocatable :: ylm (:)
      Complex (8), Allocatable :: zwork (:)
      Allocate (tp(2, ntpll))
      Allocate (rlm(lmmaxvr))
      Allocate (ylm(lmmaxvr))
      Allocate (ipiv(lmmaxvr))
      lwork = 2 * lmmaxvr
      Allocate (work(lwork))
      Allocate (zwork(lwork))

! ! allocate real SHT matrices for lmaxapw
!       If (allocated(rbshtapw)) deallocate (rbshtapw)
!       Allocate (rbshtapw(lmmaxapw, lmmaxapw))
!       If (allocated(rfshtapw)) deallocate (rfshtapw)
!       Allocate (rfshtapw(lmmaxapw, lmmaxapw))
! ! allocate real SHT matrices for lmaxvr
!       If (allocated(rbshtvr)) deallocate (rbshtvr)
!       Allocate (rbshtvr(lmmaxvr, lmmaxvr))
!       If (allocated(rfshtvr)) deallocate (rfshtvr)
!       Allocate (rfshtvr(lmmaxvr, lmmaxvr))
! ! allocate complex SHT matrices for lmaxapw
!       If (allocated(zbshtapw)) deallocate (zbshtapw)
!       Allocate (zbshtapw(lmmaxapw, lmmaxapw))
!       If (allocated(zfshtapw)) deallocate (zfshtapw)
!       Allocate (zfshtapw(lmmaxapw, lmmaxapw))

! allocate complex SHT matrices for lmaxvr
      If (allocated(zbshthf)) deallocate (zbshthf)
      Allocate (zbshthf(ntpll, lmmaxvr))
      If (allocated(zfshthf)) deallocate (zfshthf)
      Allocate (zfshthf(lmmaxvr, ntpll))

! ! generate spherical covering set for lmaxapw
!       Call sphcover (lmmaxapw, tp)

! ! generate real and complex spherical harmonics and set the backward SHT arrays
!       Do itp = 1, lmmaxapw
!          Call genrlm (input%groundstate%lmaxapw, tp(:, itp), rlm)
!          rbshtapw (itp, 1:lmmaxapw) = rlm (1:lmmaxapw)
!          Call genylm (input%groundstate%lmaxapw, tp(:, itp), ylm)
!          zbshtapw (itp, 1:lmmaxapw) = ylm (1:lmmaxapw)
!       End Do

! generate spherical covering set for lmaxvr
      Call sphcover (ntpll, tp)

#ifdef XS
      If (allocated(sphcov)) deallocate (sphcov)
      Allocate (sphcov(4, ntpll))
      If (allocated(sphcovtp)) deallocate (sphcovtp)
      Allocate (sphcovtp(2, ntpll))
      Call ld0770(sphcov (1, :),sphcov (2, :),sphcov (3, :),sphcov (4, :),ntpll)

      tp(1,:) = atan2(sqrt(sphcov(1,:)**2+sphcov(2,:)**2),sphcov(3,:))
      tp(2,:) = atan2(sphcov(2,:),sphcov(1,:))
      ! sphcovtp (:, :) = tp (:, :)
      ! sphcov (1, :) = Sin (tp(1, :)) * Cos (tp(2, :))
      ! sphcov (2, :) = Sin (tp(1, :)) * Sin (tp(2, :))
      ! sphcov (3, :) = Cos (tp(1, :))
#endif

! output data into a file
! open(1, file = 'sphcov.dat')
! do i=1,ntpll
!       ! write(1,*) sphcov (1, i),sphcov (2, i),sphcov (3, i)
!       write(1,*) Sin (tp(1, i)) * Cos (tp(2, i)),Sin (tp(1, i)) * Sin (tp(2, i)),Cos (tp(1, i))
! end do
! close(1)

! generate real and complex spherical harmonics and set the backward SHT arrays
      Do itp = 1, ntpll
         Call genylm (input%groundstate%lmaxvr, tp(:, itp), ylm)
         zbshthf (itp, 1:lmmaxvr) = ylm (1:lmmaxvr)
      End Do
! generate real and complex spherical harmonics and set the backward SHT arrays
      Do itp = 1, ntpll
         Call genylm (input%groundstate%lmaxvr, tp(:, itp), ylm)
      !    Write(*,*) "sum(ylm)", sum(ylm (1:lmmaxvr))
      !    Write(*,*) "w", sphcov(4, itp)
         zfshthf (1:lmmaxvr, itp) = fourpi * sphcov(4, itp) * conjg(ylm (1:lmmaxvr))
      End Do
End Subroutine
!EOC
