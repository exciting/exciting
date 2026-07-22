!
subroutine write_vxnl()
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use modmain
  use mod_hybrids
  use precision, only: i32, long_int
!
!   Writes the diagonal matrix elements of the non-local potential
!   into the file VXNL.OUT
!
  implicit none
  integer(long_int) :: Recl
  integer(i32) :: ik, ikfirst, iklast
  integer(i32) :: ist
  integer(i32) :: fid

!$OMP CRITICAL

  ! overwrite existing files
  open(newunit=fid, File=fname_vxnl, form='UNFORMATTED', status='REPLACE')
  close(fid)

  call inquire_large( Recl, [nkpt, nstfv], vxnl(:,:,1) )
  call open_direct_unformatted_large( fid, fname_vxnl, "write", Recl, "old" )
  do ik = 1, nkpt
      write(fid,rec=ik) nkpt, nstfv, vxnl(:,:,ik)
  end do
  close(fid)

!$OMP END CRITICAL

  return
end subroutine
