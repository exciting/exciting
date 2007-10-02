
module m_getemat
  implicit none
contains

  subroutine getemat(iq,ik,tarec,filnam,xou,xuo)
    use modmain
    use modtddft
    use modpar
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    logical :: tarec
    character(*) :: filnam
    complex(8), intent(out) :: xou(:,:,:), xuo(:,:,:)
    ! local variables
    character(*), parameter :: thisnam = 'getemat'
    integer :: recl, un, ikr, nstval_, nstcon_, nkpt_, ngq_
    real(8) :: vql_(3), vkl_(3)
    logical :: existent, opened
    ! functions
    real(8) :: r3dist
    external :: r3dist
  
    ! TO DO --------------------------------------------------
    ! selection of window of states relative to Fermi level
    ! selection of gqmax 
    ! determine record length
    ! --------------------------------------------------------

    ! check if file exists
    inquire(file=trim(filnam),exist=existent)
    if (.not.existent) then
       write(unitout,'(a)') 'Error('//thisnam//'): file does not exist: '// &
            trim(filnam)
       call terminate()
    end if

    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(nproc,nkpt,ik,ikr)
    ! I/O record length
    inquire(iolength=recl) nstval_, nstcon_, nkpt_, ngq_, vql_, vkl_, &
         xou, xuo
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted', &
         action='read', access='direct',recl=recl)

    ! read
    read(un,rec=ikr) nstval_, nstcon_, nkpt_, ngq_, vql_, vkl_, &
         xou, xuo
    close(un)

    ! check consistency
    if ( (nstval_.ne.nstval).or.(nstcon_.ne.nstcon).or.(nkpt_.ne.nkpt).or. &
         (r3dist(vql_,vql(1,iq)) > epslat).or.&
         (r3dist(vkl_,vkl(1,ik)) > epslat) ) then
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &matrix elements (current/file): '
       write(unitout,'(a,2i6)') 'nstval', nstval, nstval_
       write(unitout,'(a,2i6)') 'nstcon', nstcon, nstcon_
       write(unitout,'(a,2i6)') 'nkpt', nkpt, nkpt_
       write(unitout,'(a,3f12.6,a,3f12.6)') 'vql', vql(:,iq), ',', vql_
       write(unitout,'(a,3f12.6,a,3f12.6)') 'vkl', vkl(:,ik), ',', vkl_
       call terminate()
    end if

  end subroutine getemat

end module m_getemat
