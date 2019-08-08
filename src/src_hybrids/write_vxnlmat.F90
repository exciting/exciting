
subroutine write_vxnlmat()

  use modmain
  use mod_hybrids
  use modmpi
  use m_getunit

  implicit none
  integer(8) :: recl
  integer :: ik, ikfirst, iklast
  integer :: fid

!$OMP CRITICAL

  ! overwrite existing files
  if (rank==0) then
    call getunit(fid)
    open(fid, File=fname_vxnlmat, form='UNFORMATTED', status='REPLACE')
    close(fid)
  endif
  call barrier

  ikfirst = firstk(rank)
  iklast = lastk(rank)
  do ik = 1, nkpt
    if ((ik >= ikfirst).and.(ik <= iklast)) then ! choose a process that is going to write
      call getunit(fid)
      inquire(iolength=Recl) nkpt, nmatmax ,vnlmat(:,:,ik)
      open(fid, File=fname_vxnlmat, Action='WRITE', Form='UNFORMATTED', &
           Access='DIRECT', Status='OLD', Recl=Recl)
      write(fid, Rec=ik) nkpt, nmatmax, vnlmat(:,:,ik)
      close(fid)
    end if ! rank
    call barrier
  end do

  ! Non-local potential energy
  if (rank == 0) then
    call getunit(fid)
    open(fid, File='EXNL.OUT', Action='WRITE')
    write(fid,*) exnl
    close(fid)
  end if

!$OMP END CRITICAL
  return
end subroutine
