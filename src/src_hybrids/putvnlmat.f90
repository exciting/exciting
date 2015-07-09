!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! !REVISION HISTORY:
!   Created March 2015 (UW)
!

Subroutine putvnlmat 
      Use modmain
      use mod_hybrids
      Implicit None
  ! local variables
      Integer(8) :: recl,ik
  !$OMP CRITICAL
    Inquire (IoLength=Recl) nkpt, nmatmax ,vnlmat(:,:,1)
    Open (70, File='VNLMAT.OUT', Action='WRITE', Form='UNFORMATTED', &
   &   Access='DIRECT', status='REPLACE', Recl=Recl)
    do ik = 1, nkpt
        write(70, Rec=ik) nkpt, nmatmax ,vnlmat(:,:,ik)
    end do ! ik
    Close(70)
 !$OMP END CRITICAL
      Return
End Subroutine putvnlmat
