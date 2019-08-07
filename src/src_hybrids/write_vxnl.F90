!
subroutine write_vxnl()
  use modmain
  use mod_hybrids
  use modmpi
  use m_getunit
!
!   Writes the diagonal matrix elements of the non-local potential
!   into the file VXNL.OUT
!
  implicit none
  integer :: Recl
  integer :: ik, ikfirst, iklast
  integer :: ist
  integer :: fid
  character(80) :: fname

  ! Save < nk | \Sigma_x | nk >

!$OMP CRITICAL

  fname = 'VXNL.OUT'

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
      inquire(iolength=Recl) nkpt, nstfv, vxnl(:,:,ik)
      call getunit(fid)
      open(fid, File=fname, Action='WRITE', Form='UNFORMATTED', &
           Access='DIRECT', Status='OLD', Recl=Recl)
      write(fid,rec=ik) nkpt, nstfv, vxnl(:,:,ik)
      close(fid)
    end if ! rank
    call barrier
  end do

!$OMP END CRITICAL

  return
end subroutine
