!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genlofr
! !INTERFACE:
!
!
Subroutine genlofr
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Generates the local-orbital radial functions. This is done by integrating
!   the scalar relativistic Schr\"{o}dinger equation (or its energy deriatives)
!   at the current linearisation energies using the spherical part of the
!   effective potential. For each local-orbital, a linear combination of
!   {\tt lorbord} radial functions is constructed such that its radial
!   derivatives up to order ${\tt lorbord}-1$ are zero at the muffin-tin radius.
!   This function is normalised and the radial Hamiltonian applied to it. The
!   results are stored in the global array {\tt lofr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: np, is, ia, ias, nr, ir
      Integer :: ilo, io1, io2
      Integer :: j, l, nn, info
      Real (8) :: t1
! automatic arrays
      Real (8) :: vr (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
      Real (8) :: p0 (nrmtmax, maxlorbord), p1 (nrmtmax)
      Real (8) :: q0 (nrmtmax, maxlorbord), q1 (nrmtmax, maxlorbord)
      Real (8) :: p0s (nrmtmax), q0s (nrmtmax), q1s (nrmtmax)
      Real (8) :: hp0 (nrmtmax)
! allocatable arrays
      Integer, Allocatable :: ipiv (:)
      Real (8), Allocatable :: xa (:), ya (:)
      Real (8), Allocatable :: a (:, :), b (:), c (:)
! external functions
      Real (8) :: polynom
      External polynom
! polynomial order
      np = Max (maxlorbord+1, 4)
      Allocate (ipiv(np))
      Allocate (xa(np), ya(np), c(np))
      Allocate (a(np, np), b(np))
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            vr (1:nr) = veffmt (1, 1:nr, ias) * y00
            Do ilo = 1, nlorb (is)
               l = lorbl (ilo, is)
               Do io2 = 1, lorbord (ilo, is)
! integrate the radial Schrodinger equation
                  Call rschroddme (lorbdm(io2, ilo, is), l, 0, &
                 & lorbe(io2, ilo, ias), input%groundstate%nprad, nr, &
                 & spr(:, is), vr, nn, p0(:, io2), p1, q0(:, io2), &
                 & q1(:, io2))
! normalise radial functions
                  Do ir = 1, nr
                     fr (ir) = p0 (ir, io2) ** 2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  t1 = 1.d0 / Sqrt (Abs(gr(nr)))
                  p0 (1:nr, io2) = t1 * p0 (1:nr, io2)
                  q0 (1:nr, io2) = t1 * q0 (1:nr, io2)
                  q1 (1:nr, io2) = t1 * q1 (1:nr, io2)
! set up the matrix of radial derivatives
                  Do j = 1, np
                     ir = nr - np + j
                     xa (j) = spr (ir, is)
                     ya (j) = p0 (ir, io2) / spr (ir, is)
                  End Do
                  Do io1 = 1, lorbord (ilo, is)
                     a (io1, io2) = polynom (io1-1, np, xa, ya, c, &
                    & rmt(is))
                  End Do
               End Do
! set up the target vector
               b (:) = 0.d0
               b (lorbord(ilo, is)) = 1.d0
               Call dgesv (lorbord(ilo, is), 1, a, np, ipiv, b, np, &
              & info)
               If (info .Ne. 0) Then
                  Write (*,*)
                  Write (*, '("Error(genlofr): degenerate local-orbital&
                 & radial functions")')
                  Write (*, '(" for species ", I4)') is
                  Write (*, '(" atom ", I4)') ia
                  Write (*, '(" and local-orbital ", I4)') ilo
                  Write (*, '(" ZGESV returned INFO = ", I8)') info
                  Write (*,*)
                  Stop
               End If
! generate linear superposition of radial functions
               p0s (:) = 0.d0
               q0s (:) = 0.d0
               q1s (:) = 0.d0
               Do io1 = 1, lorbord (ilo, is)
                  t1 = b (io1)
                  p0s (1:nr) = p0s (1:nr) + t1 * p0 (1:nr, io1)
                  q0s (1:nr) = q0s (1:nr) + t1 * q0 (1:nr, io1)
                  q1s (1:nr) = q1s (1:nr) + t1 * q1 (1:nr, io1)
               End Do
! normalise radial functions
               Do ir = 1, nr
                  fr (ir) = p0s (ir) ** 2
               End Do
               Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
               t1 = 1.d0 / Sqrt (Abs(gr(nr)))
               p0s (1:nr) = t1 * p0s (1:nr)
               q0s (1:nr) = t1 * q0s (1:nr)
               q1s (1:nr) = t1 * q1s (1:nr)
! apply the scalar relativistic Hamiltonian
               Call rschrodapp (l, nr, spr(:, is), vr, p0s, q0s, q1s, &
              & hp0)
! divide by r and store in global array
               Do ir = 1, nr
                  t1 = 1.d0 / spr (ir, is)
                  lofr (ir, 1, ilo, ias) = t1 * p0s (ir)
                  lofr (ir, 2, ilo, ias) = t1 * hp0 (ir)
               End Do
            End Do
         End Do
      End Do
      Deallocate (ipiv, xa, ya, a, b, c)
      Return
End Subroutine
!EOC
