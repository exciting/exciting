!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: ggamt
! !INTERFACE:
!
!
Subroutine ggamt (is, ia, grhomt, gupmt, gdnmt, g2upmt, g2dnmt, &
& g3rhomt, g3upmt, g3dnmt)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   is      : species number (in,integer)
!   ia      : atom number (in,integer)
!   grhomt  : |grad rho| (out,real(lmmaxvr,nrmtmax))
!   gupmt   : |grad rhoup| (out,real(lmmaxvr,nrmtmax))
!   gdnmt   : |grad rhodn| (out,real(lmmaxvr,nrmtmax))
!   g2upmt  : grad^2 rhoup (out,real(lmmaxvr,nrmtmax))
!   g2dnmt  : grad^2 rhodn (out,real(lmmaxvr,nrmtmax))
!   g3rhomt : (grad rho).(grad |grad rho|) (out,real(lmmaxvr,nrmtmax))
!   g3upmt  : (grad rhoup).(grad |grad rhoup|) (out,real(lmmaxvr,nrmtmax))
!   g3dnmt  : (grad rhodn).(grad |grad rhodn|) (out,real(lmmaxvr,nrmtmax))
! !DESCRIPTION:
!   Computes $|\nabla\rho|$, $|\nabla\rho^{\uparrow}|$,
!   $|\nabla\rho^{\downarrow}|$, $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho\cdot(\nabla|\nabla\rho|)$,
!   $\nabla\rho^{\uparrow}\cdot(\nabla|\nabla\rho^{\uparrow}|)$ and
!   $\nabla\rho^{\downarrow}\cdot(\nabla|\nabla\rho^{\downarrow}|)$
!   for a muffin-tin charge density, as required by the generalised gradient
!   approximation functional for spin-polarised densities. In the case of spin
!   unpolarised calculations, $|\nabla\rho|$, $\nabla^2\rho$ and
!   $\nabla\rho\cdot(\nabla|\nabla\rho|)$ are returned in the arrays
!   {\tt gupmt}, {\tt g2upmt} and {\tt g3upmt}, respectively, while
!   {\tt grhomt}, {\tt gdnmt}, {\tt g2dnmt}, {\tt g3rhomt} and {\tt g3dnmt} are
!   not referenced. The input densities are in terms of real spherical harmonic
!   expansions but the returned functions are in spherical coordinates. See
!   routines {\tt potxc}, {\tt modxcifc}, {\tt gradrfmt}, {\tt genrlm} and
!   {\tt genshtmat}.
!
! !REVISION HISTORY:
!   Created April 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Real (8), Intent (Out) :: grhomt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: gupmt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: gdnmt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: g2upmt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: g2dnmt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: g3rhomt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: g3upmt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: g3dnmt (lmmaxvr, nrmtmax)
! local variables
      Integer :: ias, nr, ir, i, itp
! allocatable arrays
      Real (8), Allocatable :: rfmt1 (:, :)
      Real (8), Allocatable :: rftp1 (:, :, :)
      Real (8), Allocatable :: rftp2 (:, :, :)
      Real (8), Allocatable :: rftp3 (:)
      Real (8), Allocatable :: grfmt1 (:, :, :)
      Real (8), Allocatable :: grfmt2 (:, :, :)
      Allocate (rfmt1(lmmaxvr, nrmtmax))
      Allocate (rftp1(lmmaxvr, nrmtmax, 3))
      Allocate (rftp2(lmmaxvr, nrmtmax, 3))
      Allocate (rftp3(lmmaxvr))
      Allocate (grfmt1(lmmaxvr, nrmtmax, 3))
      Allocate (grfmt2(lmmaxvr, nrmtmax, 3))
      ias = idxas (ia, is)
      nr = nrmt (is)
      If (associated(input%groundstate%spin)) Then
! rhoup for spin-polarised case
         rfmt1 (:, 1:nr) = 0.5d0 * (rhomt(:, 1:nr, ias)+magmt(:, 1:nr, &
        & ias, ndmag))
      Else
! rho for spin-unpolarised case
         rfmt1 (:, 1:nr) = rhomt (:, 1:nr, ias)
      End If
! |grad rhoup|
      Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, &
     & nrmtmax, rfmt1, grfmt1)
      Do ir = 1, nr
         Do i = 1, 3
            Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, &
           & grfmt1(:, ir, i), 1, 0.d0, rftp1(:, ir, i), 1)
         End Do
         Do itp = 1, lmmaxvr
            gupmt (itp, ir) = Sqrt (rftp1(itp, ir, 1)**2+rftp1(itp, ir, &
           & 2)**2+rftp1(itp, ir, 3)**2)
         End Do
      End Do
! grad^2 rhoup
      g2upmt (:, 1:nr) = 0.d0
      Do i = 1, 3
         Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), &
        & lmmaxvr, nrmtmax, grfmt1(:, :, i), grfmt2)
         Do ir = 1, nr
            Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, &
           & grfmt2(:, ir, i), 1, 1.d0, g2upmt(:, ir), 1)
         End Do
      End Do
! (grad rhoup).(grad |grad rhoup|)
      Do ir = 1, nr
         Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, &
        & gupmt(:, ir), 1, 0.d0, rfmt1(:, ir), 1)
      End Do
      Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, &
     & nrmtmax, rfmt1, grfmt1)
      g3upmt (:, 1:nr) = 0.d0
      Do i = 1, 3
         Do ir = 1, nr
            Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, &
           & grfmt1(:, ir, i), 1, 0.d0, rftp3, 1)
            Do itp = 1, lmmaxvr
               g3upmt (itp, ir) = g3upmt (itp, ir) + rftp1 (itp, ir, i) &
              & * rftp3 (itp)
            End Do
         End Do
      End Do
      If (associated(input%groundstate%spin)) Then
! rhodn
         rfmt1 (:, 1:nr) = 0.5d0 * (rhomt(:, 1:nr, ias)-magmt(:, 1:nr, &
        & ias, ndmag))
! |grad rhodn|
         Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), &
        & lmmaxvr, nrmtmax, rfmt1, grfmt1)
         Do ir = 1, nr
            Do i = 1, 3
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, grfmt1(:, ir, i), 1, 0.d0, rftp2(:, ir, i), 1)
            End Do
            Do itp = 1, lmmaxvr
               gdnmt (itp, ir) = Sqrt (rftp2(itp, ir, 1)**2+rftp2(itp, &
              & ir, 2)**2+rftp2(itp, ir, 3)**2)
            End Do
         End Do
! grad^2 rhodn
         g2dnmt (:, 1:nr) = 0.d0
         Do i = 1, 3
            Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), &
           & lmmaxvr, nrmtmax, grfmt1(:, :, i), grfmt2)
            Do ir = 1, nr
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, grfmt2(:, ir, i), 1, 1.d0, g2dnmt(:, ir), 1)
            End Do
         End Do
! (grad rhodn).(grad |grad rhodn|)
         Do ir = 1, nr
            Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, &
           & gdnmt(:, ir), 1, 0.d0, rfmt1(:, ir), 1)
         End Do
         Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), &
        & lmmaxvr, nrmtmax, rfmt1, grfmt1)
         g3dnmt (:, 1:nr) = 0.d0
         Do i = 1, 3
            Do ir = 1, nr
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, grfmt1(:, ir, i), 1, 0.d0, rftp3, 1)
               Do itp = 1, lmmaxvr
                  g3dnmt (itp, ir) = g3dnmt (itp, ir) + rftp2 (itp, ir, &
                 & i) * rftp3 (itp)
               End Do
            End Do
         End Do
! |grad rho|
         Do ir = 1, nr
            Do itp = 1, lmmaxvr
               grhomt (itp, ir) = Sqrt ((rftp1(itp, ir, 1)+rftp2(itp, &
              & ir, 1))**2+(rftp1(itp, ir, 2)+rftp2(itp, ir, &
              & 2))**2+(rftp1(itp, ir, 3)+rftp2(itp, ir, 3))**2)
            End Do
         End Do
! (grad rho).(grad |grad rho|)
         Do ir = 1, nr
            Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, &
           & grhomt(:, ir), 1, 0.d0, rfmt1(:, ir), 1)
         End Do
         Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), &
        & lmmaxvr, nrmtmax, rfmt1, grfmt1)
         g3rhomt (:, 1:nr) = 0.d0
         Do i = 1, 3
            Do ir = 1, nr
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, grfmt1(:, ir, i), 1, 0.d0, rftp3, 1)
               Do itp = 1, lmmaxvr
                  g3rhomt (itp, ir) = g3rhomt (itp, ir) + (rftp1(itp, &
                 & ir, i)+rftp2(itp, ir, i)) * rftp3 (itp)
               End Do
            End Do
         End Do
      End If
      Deallocate (rfmt1, rftp1, rftp2, rftp3, grfmt1, grfmt2)
      Return
End Subroutine
!EOC
