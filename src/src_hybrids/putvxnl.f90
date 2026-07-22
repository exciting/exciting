
subroutine putvxnl()

  use modmain
  use mod_hybrids
  use modmpi
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use precision, only: i32, long_int, dp

  implicit none
  integer(long_int) :: recl
  integer(i32) :: ik, ikfirst, iklast, io_unit

!$OMP CRITICAL

  ikfirst = firstofset(rank, nkpt)
  iklast = lastofset(rank, nkpt)

  ! Save < m | \Sigma_x | n >
  call inquire_large( recl, [nkpt, nstfv] ,vxnl(:,:,ikfirst) )
  call open_direct_unformatted_large( io_unit, 'VXNL.OUT', "write", recl, "replace" )

  do ik = 1, nkpt
    ! check which rank should print
    if ((ik >= ikfirst).and.(ik <= iklast)) then
      write(io_unit, Rec=ik) nkpt, nstfv ,vxnl(:,:,ik)
    end if
    call barrier
  end do ! ik
  close(io_unit)

!$OMP END CRITICAL

  return
end subroutine

