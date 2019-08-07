!
!BOP
! !ROUTINE: writevnlmat
! !INTERFACE:
!
subroutine writevnlmat()
! !USES:
  use modmain
  use mod_hybrids
  use modmpi
  use m_getunit
!
! !DESCRIPTION:
!   Writes the APW matrix elements of the non-local potential
!   into the file VNLMAT.OUT
!
! !REVISION HISTORY:
!   Created March 2015 (UW)
!   Changed May 2016 (DIN)
!
!EOP
!BOC

  implicit none
  integer(8) :: recl
  integer :: ik, ikfirst, iklast
  integer :: fid
  character(80) :: fname

!$OMP CRITICAL

  fname = 'VNLMAT.OUT'

  ! overwrite existing files
  if (rank==0) then
    call getunit(fid)
    open(fid, File=fname, form='UNFORMATTED', status='REPLACE')
    close(fid)
  endif
  call barrier

  ikfirst = firstk(rank)
  iklast = lastk(rank)
  do ik = 1, nkpt
    if ((ik >= ikfirst).and.(ik <= iklast)) then ! should be the right rank ?
      inquire(iolength=Recl) nkpt, nmatmax ,vnlmat(:,:,ik)
      call getunit(fid)
      open(fid, File=fname, Action='WRITE', Form='UNFORMATTED', &
           Access='DIRECT', Status='OLD', Recl=Recl)
      write(fid, Rec=ik) nkpt, nmatmax, vnlmat(:,:,ik)
      close(fid)
    end if ! rank
    call barrier
  end do

  inquire(iolength=Recl) exnl
  call getunit(fid)
  open(fid, File=fname, Action='WRITE', Form='UNFORMATTED', &
        Access='DIRECT', Status='OLD', Recl=Recl)
  write(fid, Rec=nkpt+1) exnl
  close(fid)

!$OMP END CRITICAL
  return
end subroutine
