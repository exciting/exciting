!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putevalfv (ik, evalfv)
      Use modmain
      Use modmpi
      Implicit None
  ! arguments
      Integer, Intent (In) :: ik
      Real (8), Intent (In) :: evalfv (nstfv, nspnfv)
!
  ! local variables
      Integer :: recl, koffset
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
!
  !find the record length
      Inquire (IoLength=Recl) vkl (:, ik), nstfv, nspnfv, evalfv
!$OMP CRITICAL
      filetag = 'EVALFV'
      If (splittfile .Or. (rank .Eq. 0).or. (.not.input%sharedfs)) Then
         Open (70, File=outfilenamestring(filetag, ik), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         If (splittfile) Then
            koffset = ik - firstk (procofk(ik)) + 1
         Else
            koffset = ik
         End If
         Write (70, Rec=koffset) vkl (:, ik), nstfv, nspnfv, evalfv
         Close (70)
!
      End If
!$OMP END CRITICAL
      Return
End Subroutine putevalfv
