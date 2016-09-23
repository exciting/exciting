
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine writeprecond (ik, n, X, w)
      Use modmain
      Use modmpi
      Implicit None
      Integer, Intent (In) :: n, ik
      Complex (8), Intent (In) :: X (nmatmax, nmatmax)
      Real (8), Intent (In) :: w (nmatmax)
  !local variables
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
      Integer :: recl, koffset
      Inquire (IoLength=Recl) X, w
      filetag = "PRECONDMATRIX"
      If (splittfile .Or. (rank .Eq. 0)) Then
         Open (70, File=outfilenamestring(filetag, ik), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         If (splittfile) Then
            koffset = ik - firstk (procofk(ik)) + 1
         Else
            koffset = ik
         End If
         Write (70, Rec=koffset) X, w
         Close (70)
      Else
         Write (*,*) "Error"
         Stop
      End If
End Subroutine writeprecond
