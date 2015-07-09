
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! !REVISION HISTORY:
!   Created March 2015 (UW)
!
Subroutine getvnlmat 
      Use modmain
      Use mod_hybrids 
      Implicit None
  ! local variables
      Integer(8) :: recl,ik
      Character (256) :: fname
      logical :: exist
  !$OMP CRITICAL
    fname='VNLMAT.OUT'
    inquire(File=fname,Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(getvnlmat): File VNLMAT.OUT does not exist!'
      stop
    end if
    inquire(IoLength=Recl) nkpt, nmatmax
    open (70, File=fname, Action='READ', Form='UNFORMATTED', &
    &  Access='DIRECT', Recl=Recl)
    read(70, Rec=1) nkpt, nmatmax
    close(70)
    if (allocated(vnlmat)) deallocate(vnlmat)
    allocate(vnlmat(nmatmax,nmatmax,nkpt))
    inquire(IoLength=Recl) nkptnr, nmatmax ,vnlmat(:,:,1)
    open (70, File=fname, Action='READ', Form='UNFORMATTED', &
    &  Access='DIRECT', Recl=Recl)
    do ik = 1, nkpt
      read(70, Rec=ik) nkpt, nmatmax ,vnlmat(:,:,ik)
    end do ! ik
    close(70)
 !$OMP END CRITICAL
      Return
End Subroutine getvnlmat
