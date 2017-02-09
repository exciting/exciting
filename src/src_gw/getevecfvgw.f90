
subroutine getevecfvgw(ik,evecfv)

    use modmain
    use m_getunit
    implicit none
    ! input varibles
    integer, intent(in) :: ik
    complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
    ! local variables
    integer(8) :: recl
    integer(4) :: fid
    integer(4) :: nmatmax_, nstfv_, nspnfv_
    real(8) :: vkl_(3)
    logical :: exist
    character(256) :: filename
    
    filename = 'EVECFV_GW.OUT'
    
    inquire(File=trim(filename),Exist=exist)
    if (exist) Then
      inquire(IoLength=recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
      call getunit(fid)
      open(fid, File=trim(filename), Action='READ', Form='UNFORMATTED', &
      &    Access='DIRECT', Recl=recl)
    else
      write(*,*) 'ERROR(getevecfvgw): File ', trim(filename), ' does not exist!'
      stop
    end if
    read(fid,Rec=ik) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
    close(fid)
    
    return
end subroutine
