!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genphasedm (iq, jsym, nmax, n, phfdm, tphf)
      Use modmain
      Use modinput
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, jsym, nmax, n
      Complex (8), Intent (Out) :: phfdm (nmax, nmax)
  ! true if non-trivial phase appears at least for one (G,Gp) component
      Logical, Intent (Out) :: tphf
  ! local variables
      Real (8), Parameter :: epsortho = 1.d-12
      Real (8) :: vtl (3), t1, t2, t3
      Integer :: igq1, igq2, ivg1 (3), ivg2 (3), iv (3)
      Do igq1 = 1, n
         ivg1 (:) = ivg (:, igqig(igq1, iq))
         Do igq2 = igq1, n
        ! G-vector difference
            ivg2 (:) = ivg1 (:) - ivg (:, igqig(igq2, iq))
        ! translation vector vtl(s)
            vtl = vtlsymc (:, jsym)
            Call r3frac (input%structure%epslat, vtl, iv)
            t1 = twopi * dot_product (dble(ivg2), vtl)
            t2 = Cos (t1)
            t3 = Sin (t1)
            If (Abs(t2) .Lt. epsortho) t2 = 0.d0
            If (Abs(t3) .Lt. epsortho) t3 = 0.d0
        ! phase factor for dielectric matrix (due to translations)
            phfdm (igq1, igq2) = cmplx (t2, t3, 8)
            phfdm (igq2, igq1) = conjg (phfdm(igq1, igq2))
            If (input%xs%dbglev .Gt. 2) Then
               Write (40, '(a, i5, 2x, 2i5, 2x, i5, 2g18.10)') 'q, g, g&
              &p, jsym, phf', iq, igq1, igq2, jsym, phfdm (igq1, igq2)
            End If
        ! end loop over (G,Gp)-vectors
         End Do
      End Do
  ! occurrance of non-trivial phase for q-point
      tphf = .False.
      If (any(Abs(phfdm(1:n, 1:n)-1.d0) .Gt. input%structure%epslat)) &
     & tphf = .True.
End Subroutine genphasedm
