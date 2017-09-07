
! Copyright (C) 2005-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
! function to compose filename for parallel execution
! REMARK: never call with stringconstat allways call by reference to filetag
!
!
Character (256) Function outfilenamestring (filetag, ik)
      Use modmpi, Only: procs, lastk, firstk, procofk, splittfile
      Use modmain, Only: scrpath, task
      Use mod_misc, Only: filext
      Use modinput
      Implicit None
!external lastk,firstk
!
!character(256):: outfilenamestring
      Character (256), Intent (In) :: filetag
      Integer, Intent (In) :: ik
      Character (256) :: tmp, tmp2, krange, scrpathtmp
      krange = ''
      tmp = ''
      tmp2 = ''
      scrpathtmp = ''
      outfilenamestring = trim (scrpathtmp) // trim (filetag)// trim(filext)
#ifdef MPI
      If ((task .Eq. 0) .Or. (task .Eq. 1) .Or. (task .Eq. 2) .Or. &
      &   (task .Eq. 3) .Or. (task.eq.200) .Or. (task .Eq. 7)) Then
         If ((procs .Gt. 1) .And. splittfile) Then
            Write (tmp, '(I5)') firstk (procofk(ik))
            Write (tmp2, '(I5)') lastk (procofk(ik))
            krange = trim (adjustl(tmp)) // '-' // trim (adjustl(tmp2))
            scrpathtmp = scrpath
             outfilenamestring = trim (scrpathtmp) // trim (filetag) // trim(krange) // trim(filext)
         End If
      End If
#endif

End Function outfilenamestring
