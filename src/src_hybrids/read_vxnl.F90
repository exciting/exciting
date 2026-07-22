
subroutine read_vxnl()

    use mod_large_io, only: inquire_large, open_direct_unformatted_large
    use modmain,     only: nkpt, nstfv, wkpt
    use mod_bands,   only: nomax
    use mod_hybrids, only: vxnl, fname_vxnl
    use modmpi,      only: rank
    use precision, only: i32, long_int
    implicit none

    ! local variables
    integer(i32) :: ik, nkpt_, ib, nstfv_
    integer  :: ikfirst, iklast
    integer(i32) :: fid
    integer(long_int) :: Recl
    logical  :: exist

!$OMP CRITICAL

    inquire(File=fname_vxnl, Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(read_vxnl): File VXNL.OUT does not exist!'
      stop
    end if

    call inquire_large( Recl, [nkpt_, nstfv_] )
    call open_direct_unformatted_large( fid, fname_vxnl, "read", Recl, "old" )
    read(fid, Rec=1) nkpt_, nstfv_
    close(fid)

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

    call inquire_large( Recl, [nkpt_, nstfv_], vxnl(:,:,1) )
    call open_direct_unformatted_large( fid, fname_vxnl, "read", Recl, "old" )
    do ik = 1, nkpt
      read(fid, Rec=ik) nkpt_, nstfv_, vxnl(:,:,ik)
    end do ! ik
    close(fid)

!$OMP END CRITICAL

    return
end subroutine
