
module m_getpmat
  implicit none
contains

  subroutine getpmat(ik,vklt,tarec,filnam,pm)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: ik
    real(8), intent(in) :: vklt(:,:)
    logical :: tarec
    character(*) :: filnam
    complex(8), intent(out) :: pm(:,:,:)
    ! local variables
    character(*), parameter :: thisnam = 'getpmat'
    integer :: recl, un, ikr, nstval_, nstcon_, nkpt_
    real(8) :: vkl_(3)
    logical :: existent, opened
    ! functions
    real(8) :: r3dist
    external :: r3dist
  
    ! check if file exists
    inquire(file=trim(filnam),exist=existent)
    if (.not.existent) then
       write(unitout,'(a)') 'Error('//thisnam//'): file does not exist: '// &
            trim(filnam)
       call terminate()
    end if

    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)

    ! I/O record length
    inquire(iolength=recl) nstval_, nstcon_, nkpt_, vkl_, &
         pm
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted', &
         action='read', access='direct',recl=recl)

    ! read
    read(un,rec=ikr) nstval_, nstcon_, nkpt_, vkl_, &
         pm
    close(un)

    ! check consistency
    if ( (nstval_.ne.nstval).or.(nstcon_.ne.nstcon).or.(nkpt_.ne.nkpt).or. &
         ( r3dist(vkl_,vklt(1,ik)) > epslat ) ) then
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &matrix elements (current/file): '
       write(unitout,'(a)') 'file: '//trim(filnam)
       if (procs > 1) write(unitout,'(a,i6)') '(parallel) rank', rank
       write(unitout,'(a,i6)') 'k-point index', ik
       write(unitout,'(a,i6)') 'record position', ik
       write(unitout,'(a,2i6)') 'nstval', nstval, nstval_
       write(unitout,'(a,2i6)') 'nstcon', nstcon, nstcon_
       write(unitout,'(a,2i6)') 'nkpt', nkpt, nkpt_
       write(unitout,'(a,3f12.6,a,3f12.6)') 'vkl', vklt(:,ik), ',', vkl_
       call terminate()
    end if

  end subroutine getpmat

end module m_getpmat
