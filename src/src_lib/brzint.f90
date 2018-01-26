!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: brzint
! !INTERFACE:
!
!
Subroutine brzint (nsm, ngridk, nsk, ikmap, nw, wint, n, ld, e, f, g)
! !INPUT/OUTPUT PARAMETERS:
!   nsm    : level of smoothing for output function (in,integer)
!   ngridk : k-point grid size (in,integer(3))
!   nsk    : k-point subdivision grid size (in,integer(3))
!   ikmap  : map from grid to k-point set
!            (in,integer(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
!   nw     : number of energy divisions (in,integer)
!   wint   : energy interval (in,real(2))
!   n      : number of functions to integrate (in,integer)
!   ld     : leading dimension (in,integer)
!   e      : array of energies as a function of k-points (in,real(ld,*))
!   f      : array of weights as a function of k-points (in,real(ld,*))
!   g      : output function (out,real(nw))
! !DESCRIPTION:
!   Given energy and weight functions, $e$ and $f$, on the Brillouin zone and a
!   set of equidistant energies $\omega_i$, this routine computes the integrals
!   $$ g(\omega_i)=\frac{\Omega}{(2\pi)^3}\int_{\rm BZ} f({\bf k})
!    \delta(\omega_i-e({\bf k}))d{\bf k}, $$
!   where $\Omega$ is the unit cell volume. This is done by first interpolating
!   $e$ and $f$ on a finer $k$-point grid using the trilinear method. Then for
!   each $e({\bf k})$ on the finer grid the nearest $\omega_i$ is found and
!   $f({\bf k})$ is accumulated in $g(\omega_i)$. If the output function is
!   noisy then either {\tt nsk} should be increased or {\tt nw} decreased.
!   Alternatively, the output function can be artificially smoothed up to a
!   level given by {\tt nsm}. See routine {\tt fsmooth}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!   Improved efficiency May 2007 (Sebastian Lebegue)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: nsm
      Integer, Intent (In) :: ngridk (3)
      Integer, Intent (In) :: nsk (3)
      Integer, Intent (In) :: ikmap (0:ngridk(1)-1, 0:ngridk(2)-1, &
     & 0:ngridk(3)-1)
      Integer, Intent (In) :: nw
      Real (8), Intent (In) :: wint (2)
      Integer, Intent (In) :: n
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: e (ld,*)
      Real (8), Intent (In) :: f (ld,*)
      Real (8), Intent (Out) :: g (nw)
! local variables
      Integer :: i1, i2, i3, j1, j2, j3, k1, k2, k3, i, iw
      Integer :: i000, i001, i010, i011, i100, i101, i110, i111
      Real (8) :: p1, p2, p3, q1, q2, q3
      Real (8) :: es, fs, wd, dw, dwi, t1
! allocatable arrays
      Real (8), Allocatable :: f0 (:), f1 (:)
      Real (8), Allocatable :: e0 (:), e1 (:)
      Real (8), Allocatable :: f00 (:), f01 (:), f10 (:), f11 (:)
      Real (8), Allocatable :: e00 (:), e01 (:), e10 (:), e11 (:)
      If ((ngridk(1) .Lt. 1) .Or. (ngridk(2) .Lt. 1) .Or. (ngridk(3) &
     & .Lt. 1)) Then
         Write (*,*)
         Write (*, '("Error(brzint): ngridk < 1 : ", 3I8)') ngridk
         Write (*,*)
         Stop
      End If
      If ((nsk(1) .Lt. 1) .Or. (nsk(2) .Lt. 1) .Or. (nsk(3) .Lt. 1)) &
     & Then
         Write (*,*)
         Write (*, '("Error(brzint): nsk < 1 : ", 3I8)') nsk
         Write (*,*)
         Stop
      End If
      Allocate (f0(n), f1(n), e0(n), e1(n))
      Allocate (f00(n), f01(n), f10(n), f11(n))
      Allocate (e00(n), e01(n), e10(n), e11(n))
! length of interval
      wd = wint (2) - wint (1)
! energy step size
      dw = wd / dble (nw)
      dwi = 1.d0 / dw
      g (:) = 0.d0
#ifdef USEOMP
!!$omp parallel default( shared) private( j1, k1, j2, k2, j3, k3, i000, i001, i010, i011, i100, i101, i110, i111, i1, p1, q1, f00, f01, f10, f11, e00, e01, e10, e11, i2, p2, q2, f0, f1, e0, e1, i3, p3, q3, i, fs, es, t1, iw)
!!$omp do collapse(3)
#endif
      Do j1 = 0, ngridk (1) - 1
         Do j2 = 0, ngridk (2) - 1
            Do j3 = 0, ngridk (3) - 1
               k1 = Mod (j1+1, ngridk(1))
               k2 = Mod (j2+1, ngridk(2))
               k3 = Mod (j3+1, ngridk(3))
               i000 = ikmap (j1, j2, j3)
               i001 = ikmap (j1, j2, k3)
               i010 = ikmap (j1, k2, j3)
               i011 = ikmap (j1, k2, k3)
               i100 = ikmap (k1, j2, j3)
               i101 = ikmap (k1, j2, k3)
               i110 = ikmap (k1, k2, j3)
               i111 = ikmap (k1, k2, k3)
               Do i1 = 0, nsk (1) - 1
                  p1 = dble (i1) / dble (nsk(1))
                  q1 = 1.d0 - p1
                  f00 (:) = f (:, i000) * q1 + f (:, i100) * p1
                  f01 (:) = f (:, i001) * q1 + f (:, i101) * p1
                  f10 (:) = f (:, i010) * q1 + f (:, i110) * p1
                  f11 (:) = f (:, i011) * q1 + f (:, i111) * p1
                  e00 (:) = e (:, i000) * q1 + e (:, i100) * p1
                  e01 (:) = e (:, i001) * q1 + e (:, i101) * p1
                  e10 (:) = e (:, i010) * q1 + e (:, i110) * p1
                  e11 (:) = e (:, i011) * q1 + e (:, i111) * p1
                  Do i2 = 0, nsk (2) - 1
                     p2 = dble (i2) / dble (nsk(2))
                     q2 = 1.d0 - p2
                     f0 (:) = f00 (:) * q2 + f10 (:) * p2
                     f1 (:) = f01 (:) * q2 + f11 (:) * p2
                     e0 (:) = e00 (:) * q2 + e10 (:) * p2
                     e1 (:) = e01 (:) * q2 + e11 (:) * p2
                     Do i3 = 0, nsk (3) - 1
                        p3 = dble (i3) / dble (nsk(3))
                        q3 = 1.d0 - p3
                        Do i = 1, n
                           fs = f0 (i) * q3 + f1 (i) * p3
                           If (Abs(fs) .Gt. 1.d-10) Then
                              es = e0 (i) * q3 + e1 (i) * p3
                              t1 = (es-wint(1)) * dwi
                              iw = Nint (t1) + 1
                              If ((iw .Ge. 1) .And. (iw .Le. nw)) then
#ifdef USEOMP
!!$omp atomic update
#endif
                                 g (iw) = g (iw) + fs
#ifdef USEOMP
!!$omp end atomic
#endif
                              End If
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
#ifdef USEOMP
!!$omp end do
!!$omp end parallel
#endif
! normalise function
      t1 = dw * dble (ngridk(1)*ngridk(2)*ngridk(3)) * dble &
     & (nsk(1)*nsk(2)*nsk(3))
      t1 = 1.d0 / t1
      g (:) = t1 * g (:)
! smooth output function if required
      If (nsm .Gt. 0) Call fsmooth (nsm, nw, 1, g)
      Deallocate (f0, f1, e0, e1)
      Deallocate (f00, f01, f10, f11)
      Deallocate (e00, e01, e10, e11)
      Return
End Subroutine
!EOC
