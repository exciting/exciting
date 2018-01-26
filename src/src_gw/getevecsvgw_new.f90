
subroutine getevecsvgw_new(fname,ik,vpl,nmatmax,nstsv,nspinor,evecsv)
    use m_getunit
    implicit none
    ! arguments
    character(*) :: fname
    integer, intent(In) :: ik
    real(8), intent(In) :: vpl(3)
    integer, intent(In) :: nmatmax
    integer, intent(In) :: nstsv
    integer, intent(In) :: nspinor
    complex(8), intent(Out) :: evecsv(nmatmax,nstsv,nspinor)
    ! local variables
    logical :: exist
    integer :: recl, fid, igk, ist
    integer :: nmatmax_, nstsv_, nspinor_
    real(8) :: vpl_(3), t1
    complex(8), allocatable :: evecsv_(:,:,:)

    inquire(File=trim(fname),Exist=exist)
    if (exist) then
      call getunit(fid)
      inquire(IoLength=recl) nmatmax_, nstsv_, nspinor_, vpl_
      open(fid,File=trim(fname),Action='READ',Form='UNFORMATTED', &
      &    Access='DIRECT',Recl=recl)
      read(fid,Rec=1) nmatmax_, nstsv_, nspinor_, vpl_
      close(fid)
      if (nmatmax_<nmatmax) then
        write(*,*)
        write(*,'("ERROR(getevecsvgw_new): Wrong number of states")')
        write(*,'("  input size : ", i8)') nmatmax
        write(*,'(" ",a," : ",i8)') trim(fname), nmatmax_
        write(*,*)
        stop
      end if
      if (nstsv_<nstsv) then
        write(*,*)
        write(*,'("ERROR(getevalsvgw_new): Wrong number of states")')
        write(*,'("  input size : ", i8)') nstsv
        write(*,'(" ",a," : ",i8)') trim(fname), nstsv_
        write(*,*)
        stop
      end if
      if (nspinor_<nspinor) then
        write(*,*)
        write(*,'("ERROR(getevalsvgw_new): Wrong spinor number")')
        write(*,'("  input size : ", i8)') nspinor
        write(*,'(" ",a," : ",i8)') trim(fname), nspinor_
        write(*,*)
        stop
      end if
      allocate(evecsv_(nmatmax_,nstsv_,nspinor_))
      inquire(IoLength=recl) nmatmax_, nstsv_, nspinor_, vpl_, evecsv_
      open(fid,File=trim(fname),Action='READ',Form='UNFORMATTED', &
      &    Access='DIRECT',Recl=recl)
      read(fid,Rec=ik) nmatmax_, nstsv_, nspinor_, vpl_, evecsv_
      close(fid)
      t1 = dabs(vpl(1)-vpl_(1)) + &
      &    dabs(vpl(2)-vpl_(2)) + &
      &    dabs(vpl(3)-vpl_(3))
      if (t1 > 1.d-6) then
        write(*,*)
        write(*,'("ERROR(getevecsvgw_new): Differing vectors for k-point",i8)') ik
        write(*,'("  input vector : ", 3G18.10)') vpl
        write(*,'("  ", a, " : ", 3G18.10)') trim(fname), vpl_
        write(*,*)
        stop
      end if
      evecsv(:,:,:) = evecsv_(1:nmatmax, 1:nstsv, 1:nspinor)
    else
      write(*,*)
      write(*,*) "ERROR(getevecsvgw_new): File ", trim(fname), " does not exist!"
      write(*,*)
      stop
    end if
      
    return
end subroutine

