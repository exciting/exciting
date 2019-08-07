!
subroutine write_vxnl()
  use modmain
  use mod_hybrids
  use modmpi
!
!   Writes the diagonal matrix elements of the non-local potential
!   into the file VXNL.OUT
!
  implicit none
  integer :: Recl
  integer :: ik, ikfirst, iklast
  integer :: ist

  ! Save < nk | \Sigma_x | nk >

!$OMP CRITICAL

  ikfirst = firstk(rank)
  iklast = lastk(rank)

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

  ! open(70, File='VXNL.DAT', Action='WRITE', status='REPLACE')
  ! write(70, '(I6, " : nkpt")') nkpt
  ! write(70, '(I6, " : nstsv")') nstfv
  ! do ik = 1, nkpt
  !   ! check which rank should print
  !   if ((ik >= ikfirst).and.(ik <= iklast)) then
  !     write(70,*)
  !     write(70, '(I6, 3G18.10, " : k-point, vkl")') ik, vkl(:,ik)
  !     write(70, '(" (state, eigenvalue and occupancy below)")')
  !     do ist = 1, nstfv
  !       write(70, '(I6, 2G18.10)') ist, vxnl(ist,ist,ik)
  !     end do
  !     write(70,*)
  !   end if
  !   call barrier
  ! end do ! ik
  ! close(70)

!$OMP END CRITICAL
  return
end subroutine
