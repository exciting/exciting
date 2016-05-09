
subroutine getevalsvgw(ik,eval)
    use modmain
    use m_getunit
    implicit none
    integer, intent(in)  :: ik
    real(8), intent(out) :: eval(nstsv)
! local variables
    integer :: fid, recl
    integer :: nstsv_
    real(8) :: vkl_(3)
    character(256) :: filename
    logical :: exist
    
    filename = "EVALSV_GW.OUT"

!-------------------------------------------------------------------
!   read the KS eigenenergies (only for the irreducible k-points)
!-------------------------------------------------------------------

    inquire(File=trim(filename), Exist=exist)
    if (.not.exist) then
        write(*,*) 'ERROR(getevalsvgw): File does not exist!'
        stop
    end if

    ! find the record length
    inquire(IoLength=recl) vkl_, nstsv_
    call getunit(fid)  
    open(fid, File=trim(filename), Action='READ', &
    &  Form='UNFORMATTED', Access='DIRECT', Recl=recl)
    read(fid, Rec=1) vkl_, nstsv_
    close(fid)

    if (nstsv_ < nstsv) then
        write(*,*) 'ERROR(getevalsvgw): Different number of states!'
        write(*,*) '   Required: nstsv  =', nstsv
        write(*,*) '  Available: nstsv_ =', nstsv_
        stop
    end if

    eval(:) = 0.d0
    inquire(IoLength=recl) vkl_, nstsv_, eval
    open(fid, File=trim(filename), Action='READ', &
    &  Form='UNFORMATTED', Access='DIRECT', Recl=recl)
    read(fid, Rec=ik) vkl_, nstsv_, eval
    close(fid)
    
    return
end subroutine
