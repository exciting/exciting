!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putevecfv (ik, evecfv)
      Use modmain
      Use modmpi
      use mod_large_io, only: inquire_large, open_direct_unformatted_large
      use precision, only: i32, dp, long_int, str_256
!
      Implicit None
  ! arguments
      Integer, Intent (In) :: ik
      Complex (dp), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
  ! local variables
!
      Character (len=str_256) :: filetag
      Character (len=str_256), External :: outfilenamestring
      integer(i32) :: koffset, io_unit
      integer(long_int) ::reclength
      
      call inquire_large( reclength, vkl (:, ik), [nmatmax, nstfv, nspnfv], evecfv )

!$OMP CRITICAL
      filetag = 'EVECFV'
      If (splittfile .Or. (rank .Eq. 0).or. (.not.input%sharedfs)) Then
         call open_direct_unformatted_large( io_unit, outfilenamestring(filetag, ik), "write", reclength, "unknown" )
         If (splittfile) Then
            koffset = ik - firstofset (procofindex(ik, nkpt), nkpt) + 1
         Else
            koffset = ik
         End If
         Write (io_unit, Rec=koffset) vkl (:, ik), nmatmax, nstfv, nspnfv, &
        & evecfv
         Close (io_unit)
      End If
!$OMP END CRITICAL
End Subroutine putevecfv
