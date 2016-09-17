
subroutine putvxnl()

  use modmain
  use mod_hybrids
  use modmpi

  implicit none
  integer(8) :: recl
  integer    :: ik, ikfirst, iklast

!$OMP CRITICAL

  ikfirst = firstk(rank)
  iklast = lastk(rank)

  ! Save < m | \Sigma_x | n >
  inquire(IoLength=Recl) nkpt, nstfv ,vxnl(:,:,ikfirst)

  open(70, File='VXNL.OUT', Action='WRITE', Form='UNFORMATTED', &
  &    Access='DIRECT', status='REPLACE', Recl=Recl)

  do ik = 1, nkpt
    ! check which rank should print
    if ((ik >= ikfirst).and.(ik <= iklast)) then
      write(70, Rec=ik) nkpt, nstfv ,vxnl(:,:,ik)
    end if
    call barrier
  end do ! ik
  close(70)

!$OMP END CRITICAL

  return
end subroutine
