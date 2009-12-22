!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_dftim
!
Contains
!
!
      Subroutine dftim (iq, ik, filnam, t1, t2, t3, t4, leta)
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik
         Integer, Intent (In), Optional :: leta (4)
         Character (*), Intent (In) :: filnam
         Real (8), Intent (In) :: t1, t2, t3, t4
    ! local variables
         Integer :: un
         Call getunit (un)
         Open (un, File=trim(filnam), Action='write', Form='formatted', &
        & Position='append')
         Write (un,*)
         Write (un, '("Timings (CPU seconds) for q-point/k-point: ", 2i&
        &6)') iq, ik
         Write (un, '("  reading matrix elements	    : ", f14.4)') t1
         Write (un, '("  generating oscillators	    : ", f14.4)') t2
         Write (un, '("  updating response function	    : ", f14.4)') &
        & t3
         Write (un, '("  total			    : ", f14.4)') t4
         If (present(leta)) Then
            Write (un, '(a, i6, a, i6, a, i6, a, i6, a)') '  ETA :', &
           & leta (1), ' d ', leta (2), ' h ', leta (3), ' m ', leta &
           & (4), ' s '
         End If
         Close (un)
      End Subroutine dftim
!
End Module m_dftim
