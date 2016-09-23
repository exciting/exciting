! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !MODULE: m_xasgauntgen
! !DESCRIPTION:
!   Provides the subroutine to calculate the spin Gaunt coefficients needed in
!   the caluclation of plane wave matrix elements between core and conduction states.
!   More details are provided in Christian Vorwerk's Master thesis.
!
! !REVISION HISTORY:
!
!   Created JUNE 2015 by Christian Vorwerk
!EOP   
!BOC
Module m_xasgauntgen
      Implicit None
Contains
!
!
      Subroutine xasgauntgen (lmax2, lmax3)
         Use modmain
         Use modxs
         use modxas
         Implicit None
    ! arguments
         Integer, Intent (In) :: lmax2, lmax3
    ! local variables
         Integer :: n1, n2, l2, l3, m2, m3, lm1, lm2, lm3, &
        & lmmax2, lmmax3 
         Real (8) :: gaunt, prefac
         External :: gaunt
    ! allocate and generate complex Gaunt coefficient array
         lmmax2 = (lmax2+1) ** 2
         lmmax3 = (lmax3+1) ** 2
         prefac=1.0d0/sqrt(2.0d0)
         If (allocated(xsgntou)) deallocate (xsgntou)
         If (allocated(xsgntuo)) deallocate (xsgntuo)
         If  (allocated(xsgntoo)) deallocate (xsgntoo)
         Allocate (xsgntou(nxas,lmmax2, lmmax3))
         Allocate (xsgntuo(nxas,lmmax2, lmmax3))
         Allocate (xsgntoo(nxas,lmmax2, nxas))
         Do n1 = 1, nxas             
			Do l2 = 0, lmax2
                  Do m2 = - l2, l2
                     lm2 = idxlm (l2, m2)
                     Do l3 = 0, lmax3
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                          !resonant part
                           xsgntou (n1, lm2, lm3) = prefac*preml(n1+xasstart-1,1)*gaunt (lxas, l2, l3, &
                          & mj2ml(n1+xasstart-1,1), m2, m3)+prefac*preml(n1+xasstart-1,2)*gaunt(lxas, l2, l3,     &
                          & mj2ml(n1+xasstart-1,2), m2, m3)
                          ! anti-resonant
                           xsgntuo (n1, lm2, lm3) = prefac*preml(n1+xasstart-1,1)*gaunt (l3, l2, lxas, &
                          & m3, m2, mj2ml(n1+xasstart-1,1))+prefac*preml(n1+xasstart-1,2)*gaunt(l3, l2, lxas,     &
                          & m3, m2, mj2ml(n1+xasstart-1,2))
                          
                        End Do
                     End Do
                     Do n2 = 1, nxas
						xsgntoo(n1,lm2,n2)=preml(n1+xasstart-1,1)*preml(n2+xasstart-1,1)*gaunt (lxas, l2, lxas, &
                          & mj2ml(n1+xasstart-1,1), m2, mj2ml(n2+xasstart-1,1))+preml(n1+xasstart-1,2)*preml(n2+xasstart-1,2)*gaunt(lxas, l2, lxas,     &
                          & mj2ml(n1+xasstart-1,2), m2, mj2ml(n2+xasstart-1,2))
                     End Do
                  End Do
			End Do
         End Do
      End Subroutine xasgauntgen
!
End Module m_xasgauntgen
! EOC
