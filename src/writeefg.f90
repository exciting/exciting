!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writeefg
! !INTERFACE:
!
!
Subroutine writeefg
! !USES:
      Use modinput
      Use modmain
      Use Fox_wxml
! !DESCRIPTION:
!   Computes the electric field gradient (EFG) tensor for each atom, $\alpha$,
!   and writes it to the file {\tt EFG.OUT} along with its eigenvalues. The EFG
!   is defined by
!   $$ V^{\alpha}_{ij}\equiv\left.\frac{\partial^2 V'_{\rm C}({\bf r})}
!    {\partial{\bf r}_i\partial{\bf r}_j}\right|_{{\bf r}={\bf r}_{\alpha}}, $$
!   where $V'_{\rm C}$ is the Coulomb potential with the $l=m=0$ component
!   removed in each muffin-tin. The derivatives are computed explicitly using
!   the routine {\tt gradrfmt}.
!
! !REVISION HISTORY:
!   Created May 2004 (JKD)
!   Fixed serious problem, November 2006 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer, Parameter :: lwork = 10
      Integer :: is, ia, ias, ir, i, j, info
      Real (8) :: efg (3, 3), a (3, 3)
      Real (8) :: w (3), work (lwork)
      Type (xmlf_t), Save :: xf
! allocatable arrays
      Real (8), Allocatable :: rfmt (:, :)
      Real (8), Allocatable :: grfmt1 (:, :, :)
      Real (8), Allocatable :: grfmt2 (:, :, :)
! initialise universal variables
      Call init0
! read density and potentials from file
      Call readstate
! allocate local arrays
      Allocate (rfmt(lmmaxvr, nrmtmax))
      Allocate (grfmt1(lmmaxvr, nrmtmax, 3))
      Allocate (grfmt2(lmmaxvr, nrmtmax, 3))
      Open (50, File='EFG.OUT', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '("(electric field gradient tensor is in Cartesian coo&
     &rdinates)")')
      Call xml_OpenFile ("EFG.xml", xf, replace=.True., pretty_print=.True.)
      Call xml_newElement (xf, "EFG")
      Do is = 1, nspecies
         Call xml_newElement (xf, "species")
         Call xml_AddAttribute (xf, "n", is)
         Call xml_AddAttribute (xf, "chemicalSymbol", trim(input%structure &
                             & %speciesarray(is)%species%chemicalSymbol))
         Do ia = 1, natoms (is)
            Call xml_newElement (xf, "atom")
            Call xml_AddAttribute (xf, "n",ia)
            Call xml_newElement (xf, "EFG-tensor")
            ias = idxas (ia, is)
            Write (50,*)
            Write (50,*)
            Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') &
           & is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol), &
           & ia
! remove the l=m=0 part of the potential
            Do ir = 1, nrmt (is)
               rfmt (1, ir) = 0.d0
               rfmt (2:lmmaxvr, ir) = vclmt (2:lmmaxvr, ir, ias)
            End Do
! compute the gradient of the Coulomb potential
            Call gradrfmt (input%groundstate%lmaxvr, nrmt(is), spr(:, &
           & is), lmmaxvr, nrmtmax, rfmt, grfmt1)
            Do i = 1, 3
! compute the gradient of the gradient
               Call gradrfmt (input%groundstate%lmaxvr, nrmt(is), &
              & spr(:, is), lmmaxvr, nrmtmax, grfmt1(:, :, i), grfmt2)
               Do j = 1, 3
                  efg (i, j) = grfmt2 (1, 1, j) * y00
               End Do
            End Do
! symmetrise the EFG
            Do i = 1, 3
               Do j = i + 1, 3
                  efg (i, j) = 0.5d0 * (efg(i, j)+efg(j, i))
                  efg (j, i) = efg (i, j)
               End Do
            End Do
            Write (50,*)
            Write (50, '(" EFG tensor :")')
            Do i = 1, 3
               Write (50, '(3G18.10)') (efg(i, j), j=1, 3)
            End Do
            Write (50, '(" trace : ", G18.10)') efg (1, 1) + efg (2, 2) &
           & + efg (3, 3)
           Call xml_AddAttribute (xf, "trace", efg (1, 1) + efg (2, 2) &
           & + efg (3, 3))
            Do i=1, 3
               Call xml_newElement (xf, "line")
               Call xml_AddCharacters (xf, efg(i, :))
               Call xml_EndElement (xf, "line")
            End do
! diagonalise the EFG
            a (:, :) = efg (:, :)
            Call dsyev ('N', 'U', 3, a, 3, w, work, lwork, info)
            Write (50, '(" eigenvalues :")')
            Write (50, '(3G18.10)') w
            Call xml_newElement (xf, "eigenvalues")
            Call xml_AddCharacters (xf, w)
            Call xml_EndElement (xf, "eigenvalues")!
            Call xml_EndElement (xf, "EFG-tensor")
            Call xml_EndElement (xf, "atom")
         End Do
         Call xml_EndElement (xf, "species")
      End Do
      Call xml_EndElement (xf, "EFG")
      Call xml_Close (xf)
      Close (50)
      Write (*,*)
      Write (*, '("Info(writeefg): electric field gradient written to EFG.OUT")')
      Write (*,*)
      Deallocate (rfmt, grfmt1, grfmt2)
      Return
End Subroutine
!EOC
