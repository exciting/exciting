!
!
!
! Copyright (C) 2006-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine pmatxs2orig
      Use modmain
      Use modxs
      Use modmpi
      Use m_getunit
      Use m_getpmat
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'pmatxs2orig'
      Complex (8), Allocatable :: pm (:, :, :)
      Integer :: un, ik, recl
      If (rank == 0) Then
         Call init0
         Call init1
         Call init2
         Allocate (pm(3, nstsv, nstsv))
         Inquire (IoLength=Recl) pm
         Call getunit (un)
         Open (un, File='PMAT.OUT', Form='unformatted', Action='write', &
        & Status='replace', Access='direct', Recl=Recl)
         Do ik = 1, nkpt
            Call getpmat (ik, vkl, 1, nstsv, 1, nstsv, .True., 'PMAT_XS&
           &.OUT', pm)
            Write (un, Rec=ik) pm
         End Do
         Close (un)
         Deallocate (pm)
      End If
      Call barrier
      Write (unitout, '(a)') "Info(" // trim (thisnam) // "): conversio&
     &n of PMAT to original format finished"
End Subroutine pmatxs2orig
