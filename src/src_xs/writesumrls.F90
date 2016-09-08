!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writesumrls
      Implicit None
Contains
!
!
      Subroutine writesumrls (iq, s, fn)
         Use modmain
         Use modxs
         Use m_getunit
         Use m_writevars
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq
         Real (8), Intent (In) :: s (3)
         Character (*), Intent (In) :: fn
    ! local variables
         Character (*), Parameter :: thisnam = 'writesumrls'
         Call getunit (unit1)
         Open (unit1, File=trim(fn), Action='write')
    ! zeroth frequency moment sumrule
         Write (unit1, '(a, g18.10, a, g18.10, a)') 'zeroth frequency m&
        &oment sumrule (num. val. el.):', s (1), '(', chgval / 2.d0, ')&
        &'
    ! first frequency moment sumrule
         Write (unit1, '(a, g18.10, a, g18.10, a)') 'first frequency mo&
        &ment sumrule  (num. val. el.):', s (2), '(', chgval / 2.d0, ')&
        &'
    ! one over frequency sumrule
         Write (unit1, '(a, g18.10, a, g18.10, a)') 'pi half sumrule   &
        &                     (target):', s (3), '(', pi / 2.d0, ')'
         Close (unit1)
      End Subroutine writesumrls
!
End Module m_writesumrls
