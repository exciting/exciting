
module m_getx0
  implicit none
contains

  subroutine getx0(tp0,iq,iw,filnam,filxt,ch0,ch0wg,ch0hd)
    use modmain
    use modxs
    use m_getunit
    implicit none
    ! arguments
    logical, intent(in) :: tp0
    integer, intent(in) :: iq,iw
    character(*),intent(in) :: filnam,filxt
    complex(8), intent(out) :: ch0(:,:)
    complex(8), intent(out), optional :: ch0wg(:,:,:),ch0hd(:)
    ! local variables
    character(*), parameter :: thisnam = 'getx0'
    integer :: recl, un, ikr, ngq_
    real(8) :: vql_(3)
    logical :: existent

    ! check if file exists
    inquire(file=trim(filnam)//trim(filxt),exist=existent)
    if (.not.existent) then
       write(unitout,'(a)') 'Error('//trim(thisnam)//'): file does not exist:&
            &'// trim(filnam)//trim(filxt)
       call terminate
    end if

    ! q=0 but head or wings missing
    if (tp0.and.((.not.present(ch0wg)).or.(.not.present(ch0wg))) ) then
       write(*,*) 'Error('//trim(thisnam)//'): q=0 but head or wings missing'
       call terminate
    end if

    call getunit(un)
    ! I/O record length
    if (tp0) then
       inquire(iolength=recl) ngq(iq),vql(:,iq),ch0,ch0wg,ch0hd
       open(unit=un,file=trim(filnam)//trim(filxt),form='unformatted', &
            action='read',access='direct',recl=recl)
       read(un,rec=iw) ngq_,vql_,ch0,ch0wg,ch0hd
    else
       inquire(iolength=recl) ngq(iq),vql(:,iq),ch0
       open(unit=un,file=trim(filnam)//trim(filxt),form='unformatted', &
            action='read',access='direct',recl=recl)
       read(un,rec=iw) ngq_,vql_,ch0
    end if
    close(un)

    if ((ngq_.ne.ngq(iq)).or.(any(vql_.ne.vql(:,iq)))) then
       write(unitout,'(a)') 'Error('//trim(thisnam)//&
            &'): differring parameters for matrix elements (current/file): '
       write(unitout,'(a,2i6)') 'ngq', ngq(iq), ngq_
       write(unitout,'(a,3f12.6,a,3f12.6)') 'vql', vql(:,iq), ',', vql_
       write(unitout,'(a,i6)') 'for q-point :',iq
       write(unitout,'(a,i6)') 'for w-point :',iw
       call terminate
    end if
    
  end subroutine getx0

end module m_getx0
