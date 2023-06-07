
subroutine read_vxnl()

    use modmain,     only: nkpt, nstfv, wkpt
    use mod_hybrids, only: vxnl, fname_vxnl, nomax
    use modmpi,      only: rank
    use m_getunit
    implicit none

    ! local variables
    integer  :: ik, nkpt_, ib, nstfv_
    integer  :: ikfirst, iklast
    integer  :: fid
    integer  :: Recl
    logical  :: exist

!$OMP CRITICAL

    inquire(File=fname_vxnl, Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(read_vxnl): File VXNL.OUT does not exist!'
      stop
    end if

    inquire(IoLength=recl) nkpt_, nstfv_
    call getunit(fid)
    open(fid, File=fname_vxnl, Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    read(fid, Rec=1) nkpt_, nstfv_
    close(fid)

    ! consistency check
    if (nkpt_ /= nkpt) then
      write(*,*) "ERROR(read_vxnl): Inconsistent number of k-points"
      write(*,*) "nkpt  = ", nkpt
      write(*,*) "nkpt_ = ", nkpt_
      !stop
    end if
    if (nstfv_ /= nstfv) then
      write(*,*) "ERROR(read_vxnl): Inconsistent number of states"
      write(*,*) "nstfv  = ", nstfv
      write(*,*) "nstfv_ = ", nstfv_
      !stop
    end if

    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv_,nstfv_,nkpt))
    vxnl(:,:,:) = 0d0
    inquire(IoLength=recl) nkpt_, nstfv_, vxnl(:,:,1)
    call getunit(fid)
    open(fid, File=fname_vxnl, Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=Recl)
    do ik = 1, nkpt
      read(fid, Rec=ik) nkpt_, nstfv_, vxnl(:,:,ik)
    end do ! ik
    close(fid)

!$OMP END CRITICAL

    return
end subroutine
