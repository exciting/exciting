
subroutine read_vxnl()

    use modmain,     only: nkpt, nstfv
    use mod_hybrids, only: vxnl
    implicit none

    ! local variables
    integer       :: ik, nkpt_, nstfv_
    integer       :: ikfirst, iklast
    integer       :: Recl
    character(40) :: fname
    logical       :: exist

!$OMP CRITICAL

    fname = 'VXNL.OUT'

    inquire(File=trim(fname), Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(read_vxnl): File VXNL.OUT does not exist!'
      stop
    end if

    inquire(IoLength=recl) nkpt_, nstfv_
    open(70, File=fname, Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    read(70, Rec=1) nkpt_, nstfv_
    close(70)

    ! consistency check
    if (nkpt_ /= nkpt) then
      write(*,*) "ERROR(read_vxnl): Inconsistent number of k-points"
      write(*,*) "nkpt  = ", nkpt
      write(*,*) "nkpt_ = ", nkpt_
      stop
    end if
    if (nstfv_ /= nstfv) then
      write(*,*) "ERROR(read_vxnl): Inconsistent number of states"
      write(*,*) "nstfv  = ", nstfv
      write(*,*) "nstfv_ = ", nstfv_
      stop
    end if

    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,nkpt))

    inquire(IoLength=recl) nkpt_, nstfv_, vxnl(:,:,1)
    open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=Recl)
    do ik = 1, nkpt
      read(70, Rec=ik) nkpt_, nstfv_, vxnl(:,:,ik)
    end do ! ik
    close(70)

!$OMP END CRITICAL

    return
end subroutine