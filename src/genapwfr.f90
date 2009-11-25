!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genapwfr
! !INTERFACE:
!
!
Subroutine genapwfr
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Generates the APW radial functions. This is done by integrating the scalar
!   relativistic Schr\"{o}dinger equation (or its energy deriatives) at the
!   current linearisation energies using the spherical part of the effective
!   potential. The number of radial functions at each $l$-value is given by the
!   variable {\tt apword} (at the muffin-tin boundary, the APW functions have
!   continuous derivatives up to order ${\tt apword}-1$). Within each $l$, these
!   functions are orthonormalised with the Gram-Schmidt method. The radial
!   Hamiltonian is applied to the orthonormalised functions and the results are
!   stored in the global array {\tt apwfr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, nr, ir
      Integer :: nn, l, io1, io2
      Real (8) :: t1
! automatic arrays
      Real (8) :: vr (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
      Real (8) :: p0 (nrmtmax, apwordmax), p1 (nrmtmax), p1s &
     & (apwordmax)
      Real (8) :: q0 (nrmtmax, apwordmax), q1 (nrmtmax, apwordmax)
      Real (8) :: hp0 (nrmtmax)
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            vr (1:nr) = veffmt (1, 1:nr, ias) * y00
            Do l = 0, input%groundstate%lmaxapw
               Do io1 = 1, apword (l, is)
! integrate the radial Schrodinger equation
                  Call rschroddme (apwdm(io1, l, is), l, 0, apwe(io1, &
                 & l, ias), input%groundstate%nprad, nr, spr(:, is), &
                 & vr, nn, p0(:, io1), p1, q0(:, io1), q1(:, io1))
! normalise radial functions
                  Do ir = 1, nr
                     fr (ir) = p0 (ir, io1) ** 2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  t1 = 1.d0 / Sqrt (Abs(gr(nr)))
                  p0 (1:nr, io1) = t1 * p0 (1:nr, io1)
                  p1s (io1) = t1 * p1 (nr)
                  q0 (1:nr, io1) = t1 * q0 (1:nr, io1)
                  q1 (1:nr, io1) = t1 * q1 (1:nr, io1)
! subtract linear combination of previous vectors
                  Do io2 = 1, io1 - 1
                     Do ir = 1, nr
                        fr (ir) = p0 (ir, io1) * p0 (ir, io2)
                     End Do
                     Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                     t1 = gr (nr)
                     p0 (1:nr, io1) = p0 (1:nr, io1) - t1 * p0 (1:nr, &
                    & io2)
                     p1s (io1) = p1s (io1) - t1 * p1s (io2)
                     q0 (1:nr, io1) = q0 (1:nr, io1) - t1 * q0 (1:nr, &
                    & io2)
                     q1 (1:nr, io1) = q1 (1:nr, io1) - t1 * q1 (1:nr, &
                    & io2)
                  End Do
! normalise radial functions
                  Do ir = 1, nr
                     fr (ir) = p0 (ir, io1) ** 2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  t1 = Abs (gr(nr))
                  If (t1 .Lt. 1.d-20) Then
                     Write (*,*)
                     Write (*, '("Error(genapwfr): degenerate APW radia&
                    &l functions")')
                     Write (*, '(" for species ", I4)') is
                     Write (*, '(" atom ", I4)') ia
                     Write (*, '(" angular momentum ", I4)') l
                     Write (*, '(" and order ", I4)') io1
                     Write (*,*)
                     Stop
                  End If
                  t1 = 1.d0 / Sqrt (t1)
                  p0 (1:nr, io1) = t1 * p0 (1:nr, io1)
                  p1s (io1) = t1 * p1s (io1)
                  q0 (1:nr, io1) = t1 * q0 (1:nr, io1)
                  q1 (1:nr, io1) = t1 * q1 (1:nr, io1)
! apply the Hamiltonian
                  Call rschrodapp (l, nr, spr(:, is), vr, p0(:, io1), &
                 & q0(:, io1), q1(:, io1), hp0)
! divide by r and store in global array
                  Do ir = 1, nr
                     t1 = 1.d0 / spr (ir, is)
                     apwfr (ir, 1, io1, l, ias) = t1 * p0 (ir, io1)
                     apwfr (ir, 2, io1, l, ias) = t1 * hp0 (ir)
                  End Do
! derivative at the muffin-tin surface
                  apwdfr (io1, l, ias) = (p1s(io1)-p0(nr, io1)*t1) * t1
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
