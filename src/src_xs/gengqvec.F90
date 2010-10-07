!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gengqvec
! !INTERFACE:
!
!
Subroutine gengqvec (iq, vpl, vpc, ngp, igpig, vgpl, vgpc, gpc, tpgpc)
! !USES:
      Use modinput
      Use modmain
      Use modxs
! !INPUT/OUTPUT PARAMETERS:
!   vpl   : p-point vector in lattice coordinates (in,real(3))
!   vpc   : p-point vector in Cartesian coordinates (in,real(3))
!   ngp   : number of G+p-vectors returned (out,integer)
!   igpig : index from G+p-vectors to G-vectors (out,integer(ngkmax))
!   vgpl  : G+p-vectors in lattice coordinates (out,real(3,ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (out,real(3,ngkmax))
!   gpc   : length of G+p-vectors (out,real(ngkmax))
!   tpgpc : (theta, phi) coordinates of G+p-vectors (out,real(2,ngkmax))
! !DESCRIPTION:
!   Generates a set of ${\bf G+p}$-vectors for the input ${\bf p}$-point with
!   length less than {\tt gkmax}. These are used as the plane waves in the APW
!   functions. Also computes the spherical coordinates of each vector.
!   Based on {\tt gengpvec}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
      Real (8), Intent (In) :: vpl (3)
      Real (8), Intent (In) :: vpc (3)
      Integer, Intent (Out) :: ngp
      Integer, Intent (Out) :: igpig (ngqmax)
      Real (8), Intent (Out) :: vgpl (3, ngqmax)
      Real (8), Intent (Out) :: vgpc (3, ngqmax)
      Real (8), Intent (Out) :: gpc (ngqmax)
      Real (8), Intent (Out) :: tpgpc (2, ngqmax)
  ! local variables
      Integer :: ig, igp
      Real (8) :: v (3), t1, t2
!
      Integer :: isym, lspl, igpt, ivlt (3)
      Real (8) :: vl (3), vc (3), vlt (3), vct (3), vctl (3), s (3, 3), &
     & c (3, 3)
!
      If (input%xs%gqmax .Lt. input%structure%epslat) Then
         igp = 1
         igpig (igp) = igp
         vgpl (:, igp) = vpl (:)
         vgpc (:, igp) = vpc (:)
         Call sphcrd (vgpc(1, igp), gpc(igp), tpgpc(1, igp))
         ivgigq (0, 0, 0, iq) = igp
         ngp = 1
         Return
      End If
      t1 = input%xs%gqmax ** 2
      ivgigq (:, :, :, iq) = 0
      igp = 0
      Do ig = 1, ngvec
         v (:) = vgc (:, ig) + vpc (:)
         ! cutoff type for G-vectors
         if (tgqmaxg) then
           t2 = vgc(1,ig) ** 2 + vgc(2,ig) ** 2 + vgc(3,ig) ** 2
         else
           t2 = v (1) ** 2 + v (2) ** 2 + v (3) ** 2
         end if
         If (t2 .Lt. t1) Then
            igp = igp + 1
            If (igp .Gt. ngqmax) Then
               Write (*,*)
               Write (*, '("Error(gengqvec): number of G+p-vectors exce&
              &eds ngqmax")')
               Write (*,*)
               Stop
            End If
        ! index to G-vector
            igpig (igp) = ig
        ! G+p-vector in lattice coordinates
            vgpl (:, igp) = dble (ivg(:, ig)) + vpl (:)
        ! G+p-vector in Cartesian coordinates
            vgpc (:, igp) = v (:)
        ! G+p-vector length and (theta, phi) coordinates
            Call sphcrd (vgpc(1, igp), gpc(igp), tpgpc(1, igp))
        ! map from grid to G+p-vector
            ivgigq (ivg(1, ig), ivg(2, ig), ivg(3, ig), iq) = igp
         End If
      End Do
      ngp = igp
      If (input%xs%dbglev .Gt. 1) Then
         Write (*, '(a)') 'Debug(gengqvec): igp,isym,lspl,vl,vlt'
         Do igp = 1, ngp
            vl (:) = dble (ivg(:, igpig(igp)))
            vc = matmul (bvec, vl)
            Do isym = 1, nsymcrys
               lspl = lsplsymc (isym)
               c (:, :) = symlatc (:, :, lspl)
               s (:, :) = dble (symlat(:, :, lspl))
               vlt = matmul (vl, s)
               ivlt = Nint (vlt)
               vct = matmul (vc, c)
               vctl = matmul (binv, vct)
               igpt = ivgigq (ivlt(1), ivlt(2), ivlt(3), iq)
               Write (*, '(3i6,5x,3i5,3x,3i5)') igp, isym, lspl, Nint &
              & (vl), Nint (vlt)
            End Do
         End Do
         Write (*,*)
      End If
      Return
End Subroutine gengqvec
!EOC
