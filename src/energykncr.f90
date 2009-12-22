!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine energykncr
      Use modmain
      Implicit None
      Integer :: is, ia, ias, ist, ir
! allocatable local arrays
      Real (8), Allocatable :: rfmt (:, :)
! external functions
      Real (8) :: rfmtinp
      External rfmtinp
! allocate local arrays
      Allocate (rfmt(lmmaxvr, nrmtmax))
! calculate the kinetic energy for core states
      engykncr = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! sum of core eigenvalues
            Do ist = 1, spnst (is)
               If (spcore(ist, is)) engykncr = engykncr + spocc (ist, &
              & is) * evalcr (ist, ias)
            End Do
! core density
            Do ir = 1, nrmt (is)
               rfmt (1, ir) = rhocr (ir, ias) / y00
            End Do
            engykncr = engykncr - rfmtinp (1, 0, nrmt(is), spr(:, is), &
           & lmmaxvr, rfmt, veffmt(:, :, ias))
         End Do
      End Do
      Deallocate (rfmt)
      Return
End Subroutine
