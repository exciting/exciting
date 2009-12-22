!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rfarray (lmax, ld, rfmt, rfir, np, vpl, fp)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: rfmt (ld, nrmtmax, natmtot)
      Real (8), Intent (In) :: rfir (ngrtot)
      Integer, Intent (In) :: np
      Real (8), Intent (In) :: vpl (3, np)
      Real (8), Intent (Out) :: fp (np)
! local variables
      Integer :: ia, is, ias, ip, iv (3)
      Integer :: i1, i2, i3, ir0, ir, np2
      Integer :: l, m, lm, ig, ifg, i, j
      Real (8) :: rmt2, r, tp (2), sum, t1, t2
      Real (8) :: v1 (3), v2 (3), v3 (3), v4 (3), v5 (3)
! allocatable arrays
      Real (8), Allocatable :: ya (:), c (:)
      Real (8), Allocatable :: rlm (:)
      Complex (8), Allocatable :: zfft (:)
! external functions
      Real (8) :: polynom
      External polynom
      Allocate (ya(input%groundstate%nprad), &
     & c(input%groundstate%nprad))
      Allocate (rlm((lmax+1)**2))
      Allocate (zfft(ngrtot))
      np2 = input%groundstate%nprad / 2
! Fourier transform rfir to G-space
      zfft (:) = rfir (:)
      Call zfftifc (3, ngrid,-1, zfft)
! begin loop over all points
      Do ip = 1, np
         v2 (:) = vpl (:, ip)
         Call r3frac (input%structure%epslat, v2, iv)
! convert point to Cartesian coordinates
         Call r3mv (input%structure%crystal%basevect, v2, v1)
! check if point is in a muffin-tin
         Do is = 1, nspecies
            rmt2 = rmt (is) ** 2
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               v2 (:) = v1 (:) - atposc (:, ia, is)
               Do i1 = - 1, 1
                  v3 (:) = v2 (:) + dble (i1) * &
                 & input%structure%crystal%basevect(:, 1)
                  Do i2 = - 1, 1
                     v4 (:) = v3 (:) + dble (i2) * &
                    & input%structure%crystal%basevect(:, 2)
                     Do i3 = - 1, 1
                        v5 (:) = v4 (:) + dble (i3) * &
                       & input%structure%crystal%basevect(:, 3)
                        t1 = v5 (1) ** 2 + v5 (2) ** 2 + v5 (3) ** 2
                        If (t1 .Lt. rmt2) Then
                           Call sphcrd (v5, r, tp)
                           Call genrlm (lmax, tp, rlm)
                           Do ir = 1, nrmt (is)
                              If (spr(ir, is) .Ge. r) Then
                                 If (ir .Le. np2) Then
                                    ir0 = 1
                                 Else If (ir .Gt. nrmt(is)-np2) Then
                                    ir0 = nrmt (is) - &
                                   & input%groundstate%nprad + 1
                                 Else
                                    ir0 = ir - np2
                                 End If
                                 r = Max (r, spr(1, is))
                                 sum = 0.d0
                                 Do l = 0, lmax
                                    Do m = - l, l
                                       lm = idxlm (l, m)
                                       Do j = 1, &
                                      & input%groundstate%nprad
                                          i = ir0 + j - 1
                                          ya (j) = rfmt (lm, i, ias)
                                       End Do
                                       t2 = polynom (0, &
                                      & input%groundstate%nprad, &
                                      & spr(ir0, is), ya, c, r)
                                       sum = sum + t2 * rlm (lm)
                                    End Do
                                 End Do
                                 Go To 10
                              End If
                           End Do
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
! otherwise use interstitial function
         sum = 0.d0
         Do ig = 1, ngvec
            ifg = igfft (ig)
            t1 = vgc (1, ig) * v1 (1) + vgc (2, ig) * v1 (2) + vgc (3, &
           & ig) * v1 (3)
            sum = sum + dble (zfft(ifg)*cmplx(Cos(t1), Sin(t1), 8))
         End Do
10       Continue
         fp (ip) = sum
      End Do
      Deallocate (rlm, zfft, ya, c)
      Return
End Subroutine
