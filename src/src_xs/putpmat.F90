
module m_putpmat
  implicit none
contains

  subroutine putpmat(ik,tarec,filnam,pm)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: ik
    ! true if absolut record position is ik
    logical, intent(in) :: tarec
    character(*), intent(in) :: filnam
    complex(8), intent(in) :: pm(:,:,:)
    ! local variables
    integer :: un, recl, ikr
    
    ! record position for k-point
    ikr=ik
    ! record position is not absolute k-point index
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)

    ! I/O record length
    inquire(iolength=recl) nstval, nstcon, nkpt, &
         vkl(:,ik), pm
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted',action='write', &
         access='direct',recl=recl)
    write(un,rec=ikr) nstval, nstcon, nkpt, vkl(:,ik),pm
    close(un)

  end subroutine putpmat

end module m_putpmat
