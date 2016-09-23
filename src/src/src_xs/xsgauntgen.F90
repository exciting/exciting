!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_xsgauntgen
      Implicit None
Contains
!
!
      Subroutine xsgauntgen (lmax1, lmax2, lmax3)
         Use modmain
         Use modxs
         Implicit None
    ! arguments
         Integer, Intent (In) :: lmax1, lmax2, lmax3
    ! local variables
         Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3, lmmax1, &
        & lmmax2, lmmax3
         Real (8) :: gaunt
         External :: gaunt
    ! allocate and generate complex Gaunt coefficient array
         lmmax1 = (lmax1+1) ** 2
         lmmax2 = (lmax2+1) ** 2
         lmmax3 = (lmax3+1) ** 2
         If (allocated(xsgnt)) deallocate (xsgnt)
         Allocate (xsgnt(lmmax1, lmmax2, lmmax3))
         Do l1 = 0, lmax1
            Do m1 = - l1, l1
               lm1 = idxlm (l1, m1)
               Do l2 = 0, lmax2
                  Do m2 = - l2, l2
                     lm2 = idxlm (l2, m2)
                     Do l3 = 0, lmax3
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           xsgnt (lm1, lm2, lm3) = gaunt (l1, l2, l3, &
                          & m1, m2, m3)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Subroutine xsgauntgen
!
End Module m_xsgauntgen
