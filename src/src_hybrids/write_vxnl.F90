!
subroutine write_vxnl()
  use modmain
  use mod_hybrids
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

!$OMP CRITICAL

  ! overwrite existing files
  call getunit(fid)
  open(fid, File=fname_vxnl, form='UNFORMATTED', status='REPLACE')
  close(fid)

  inquire(iolength=Recl) nkpt, nstfv, vxnl(:,:,ik)
  call getunit(fid)
  open(fid, File=fname_vxnl, Action='WRITE', Form='UNFORMATTED', &
       Access='DIRECT', Status='OLD', Recl=Recl)
  do ik = 1, nkpt
      write(fid,rec=ik) nkpt, nstfv, vxnl(:,:,ik)
  end do
  close(fid)

!$OMP END CRITICAL

  return
end subroutine
