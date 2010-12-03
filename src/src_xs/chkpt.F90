!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine chkpt (ncpt, cptv, mesg)
      Use modxs
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: ncpt, cptv (ncpt)
      Character (*), Intent (In) :: mesg
#ifdef RESUME
  ! local variables
      Integer :: un
      Character (256) :: str
      Write (str,*) ncpt
      Call getunit (un)
      Open (un, File=trim(fnresume), Form='formatted', Action='write', &
     & Status='replace')
      Write (un, '(i8, " : length of checkpoint vector")') ncpt
      Write (un, '('//trim(adjustl(str))//'i8, " : checkpoint vector")') cptv (:)
      Write (un, '(" (", a, ")")') trim (mesg)
      Close (un)
#endif
End Subroutine chkpt
