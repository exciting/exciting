!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putbsediag (fname)
      Use modxs
      Use m_getunit
      use mod_constants, only: h2ev
      Implicit None
  ! arguments
      Character (*), Intent (In) :: fname
  ! local variables
      Integer :: un
      Call getunit (un)
      Open (un, File=trim(fname), Action='write', Form='formatted', &
     & Status='replace')
      Write (un, '(2g18.10, " : BSE kernel diagonal mean value")') bsed
      Write (un, '(2g18.10, " : BSE kernel diagonal lower limit")') &
     & bsedl
      Write (un, '(2g18.10, " : BSE kernel diagonal upper limit")') &
     & bsedu
      Write (un, '(2g18.10, " : BSE kernel diagonal window size")') &
     & bsedd
      Write (un,*)
      Write (un, '(2g18.10, " : BSE kernel diagonal mean value (eV)")') &
     & bsed * h2ev
      Write (un, '(2g18.10, " : BSE kernel diagonal lower limit (eV)")') bsedl * h2ev
      Write (un, '(2g18.10, " : BSE kernel diagonal upper limit (eV)")') bsedu * h2ev
      Write (un, '(2g18.10, " : BSE kernel diagonal window size (eV)")') bsedd * h2ev
      Write (un,*)
      Write (un, '(g18.10, " : BSE kernel diagonal deviation (%)")') &
     & dble (bsedd) / dble (bsed) * 100.d0
      Close (un)
End Subroutine putbsediag
