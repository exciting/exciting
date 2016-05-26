
subroutine getvnlmat()

    use modmain
    use mod_hybrids
    use modmpi
    implicit none
  
    ! local variables
    integer       :: ik, nkpt_, nmatmax_
    integer       :: ikfirst, iklast
    integer(8)    :: recl
    character(40) :: fname
    logical       :: exist

!$OMP CRITICAL

    fname = 'VNLMAT.OUT'

    inquire(File=trim(fname), Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(getvnlmat): File VNLMAT.OUT does not exist!'
      stop
    end if

    inquire(IoLength=recl) nkpt_, nmatmax_
    open(70, File=fname, Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    read(70, Rec=1) nkpt_, nmatmax_
    close(70)

    ! consistency check
    if (nkpt_ /= nkpt) then
      write(*,*) "ERROR(getvnlmat): Inconsistent number of k-points"
      write(*,*) "nkpt  = ", nkpt
      write(*,*) "nkpt_ = ", nkpt_
      stop
    end if
    if (nmatmax_ /= nmatmax) then
      write(*,*) "ERROR(getvnlmat): Inconsistent number of states"
      write(*,*) "nmatmax  = ", nmatmax
      write(*,*) "nmatmax_ = ", nmatmax_
      stop
    end if

    ikfirst = firstk(rank)
    iklast = lastk(rank)
 
     if (allocated(vnlmat)) deallocate(vnlmat)
    allocate(vnlmat(nmatmax,nmatmax,ikfirst:iklast))

    inquire(IoLength=recl) nkpt_, nmatmax_, vnlmat(:,:,ikfirst)
    open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
    &    Access='DIRECT', Recl=recl)
    do ik = 1, nkpt
      if ((ik >= ikfirst).and.(ik <= iklast)) then
        read(70, Rec=ik) nkpt_, nmatmax_, vnlmat(:,:,ik)
      end if
      call barrier
    end do ! ik
    close(70)

!$OMP END CRITICAL

    return
end subroutine
