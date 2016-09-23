!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putematrad (iqr, iq)
      Use modmain
      Use modxs
      Use m_genfilname
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: iqr, iq
  ! local variables
      Character (256) :: fname
      Integer :: un
  ! calculate radial integrals
      Call genfilname (basename='EMATRAD', iq=iqr, filnam=fname)
      Call ematrad (iq)
      Call getunit (un)
      Open (un, File=trim(fname), Form='unformatted', Action='write', &
     & Status='replace')
      Write (un) riaa, riloa, rilolo
      Close (un)
End Subroutine putematrad
