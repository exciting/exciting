
module m_puttetcw
  implicit none
contains

  subroutine puttetcw(iq,ik,iv,ic,filnam,cw,cwa,cwsurf)
    use modmain
    use modtddft
    use modpar
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,iv,ic
    character(*), intent(in) :: filnam
    real(8), intent(in) :: cw(:),cwa(:),cwsurf(:)
    ! local variables
    integer :: un, recl, irec

    ! record position
    irec=(ik-1)*nstval*nstcon + (iv-1)*nstcon + ic

    ! I/O record length
    inquire(iolength=recl) cw,cwa,cwsurf
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted', &
         action='write',access='direct',recl=recl)
    write(un,rec=irec) cw,cwa,cwsurf
    close(un)

  end subroutine puttetcw

end module m_puttetcw
