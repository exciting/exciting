
module m_putdevalsv
  implicit none
contains

  subroutine putdevalsv(iq,ik,tarec,filnam,eou,euo)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    logical :: tarec
    character(*) :: filnam
    real(8), intent(in) :: eou(:,:), euo(:,:)
    ! local variables
    integer :: un, ikr, recl
    
    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)

    ! I/O record length
    inquire(iolength=recl) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), &
         vkl(:,ik), eou, euo
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted', &
         action='write',access='direct',recl=recl)
    write(un,rec=ikr) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), vkl(:,ik), &
         eou, euo
    close(un)

  end subroutine putdevalsv

end module m_putdevalsv
