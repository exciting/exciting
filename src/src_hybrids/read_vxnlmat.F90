
subroutine read_vxnlmat()

    use modmain
    use mod_hybrids
    use modmpi
    use m_getunit

    ! local variables
    implicit none
    integer       :: ik, nkpt_, nmatmax_
    integer       :: ikfirst, iklast
    integer       :: fid
    integer(8)    :: recl
    logical       :: exist

!$OMP CRITICAL

    inquire(File=fname_vxnlmat, Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(read_vxnlmat): File VNLMAT.OUT does not exist!'
      stop
    end if

    inquire(IOLength=recl) nkpt_, nmatmax_
    call getunit(fid)
    open(fid, File=fname_vxnlmat, Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    read(fid, Rec=1) nkpt_, nmatmax_
    close(fid)

    ! consistency check
    if (nkpt_ /= nkpt) then
      write(*,*) "ERROR(read_vxnlmat): Inconsistent number of k-points"
      write(*,*) "nkpt  = ", nkpt
      write(*,*) "nkpt_ = ", nkpt_
      stop
    end if
    if (nmatmax_ /= nmatmax) then
      write(*,*) "ERROR(read_vxnlmat): Inconsistent number of states"
      write(*,*) "nmatmax  = ", nmatmax
      write(*,*) "nmatmax_ = ", nmatmax_
      stop
    end if

    ikfirst = firstk(rank)
    iklast = lastk(rank)

    if (allocated(vnlmat)) deallocate(vnlmat)
    allocate(vnlmat(nmatmax,nmatmax,ikfirst:iklast))

    call getunit(fid)
    inquire(IOLength=recl) nkpt_, nmatmax_, vnlmat(:,:,ikfirst)
    open(fid, File=trim(fname_vxnlmat), Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    do ik = ikfirst, iklast
      read(fid, Rec=ik) nkpt_, nmatmax_, vnlmat(:,:,ik)
    end do ! ik
    close(fid)

    ! Read the non-local potential  energy
    open(fid, File='EXNL.OUT', Action='READ', Status='OLD')
    read(fid,*) exnl
    close(fid)

!$OMP END CRITICAL

    return
end subroutine
