
subroutine getvxnl()

    use modmain
    use mod_hybrids 
    use modmpi
    implicit none
  
    ! local variables
    integer       :: ik, nkpt_, nstfv_
    integer       :: ikfirst, iklast
    integer(8)    :: recl
    character(40) :: fname
    logical       :: exist

!$OMP CRITICAL

    fname = 'VXNL.OUT'

    inquire(File=trim(fname), Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(getvxnl): File VXNL.OUT does not exist!'
      stop
    end if

    inquire(IoLength=recl) nkpt_, nstfv_
    open(70, File=fname, Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    read(70, Rec=1) nkpt_, nstfv_
    close(70)

    ! consistency check
    if (nkpt_ /= nkpt) then
      write(*,*) "ERROR(getvxnl): Inconsistent number of k-points"
      write(*,*) "nkpt  = ", nkpt
      write(*,*) "nkpt_ = ", nkpt_
      stop
    end if
    if (nstfv_ /= nstfv) then
      write(*,*) "ERROR(getvxnl): Inconsistent number of states"
      write(*,*) "nstfv  = ", nstfv
      write(*,*) "nstfv_ = ", nstfv_
      stop
    end if

    ikfirst = firstk(rank)
    iklast = lastk(rank)
 
    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,ikfirst:iklast))

    inquire(IoLength=recl) nkpt_, nstfv_, vxnl(:,:,ikfirst)

    open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    do ik = 1, nkpt
      if ((ik >= ikfirst).and.(ik <= iklast)) then
        read(70, Rec=ik) nkpt_, nstfv_, vxnl(:,:,ik)
      end if
    end do ! ik
    close(70)

!$OMP END CRITICAL

    return
end subroutine
