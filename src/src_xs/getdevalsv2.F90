
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getdevalsv2
  implicit none
contains

  subroutine getdevalsv2(iq,ik,tarec,filnam,e1,o1,e2,o2)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    logical :: tarec
    character(*) :: filnam
    real(8), intent(out) :: e1(:,:),o1(:,:)
    real(8), optional, intent(out) :: e2(:,:),o2(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'getdevalsv2'
    integer :: recl, un, ikr, nstval_, nstcon_, nkpt_, ngq_
    real(8) :: vql_(3), vkl_(3)
    logical :: existent
    ! functions
    real(8) :: r3dist
    external :: r3dist
    ! check if file exists
    inquire(file=trim(filnam),exist=existent)
    if (.not.existent) then
       write(unitout,'(a)') 'Error('//trim(thisnam)//'): file does not exist&
            &: '// trim(filnam)
       call terminate
    end if
    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)
    call getunit(un)
    if (present(e2).and.(.not.present(o2))) then
       write(*,*)
       write(*,'("Error(getdevalsv2): either type 1 or type 2 arguments &
            required for eigenvalue and occupancy differences - refer to &
            developers manual")')
       write(*,*)
    end if
    if (present(e2)) then
       ! I/O record length
       inquire(iolength=recl) nstval_,nstcon_,nkpt_,ngq_,vql_,vkl_, &
            e1,o1,e2,o2
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='read', access='direct',recl=recl)
       read(un,rec=ikr) nstval_,nstcon_,nkpt_,ngq_,vql_,vkl_, &
            e1,o1,e2,o2
    else
       ! I/O record length
       inquire(iolength=recl) nstval_,nstcon_,nkpt_,ngq_,vql_,vkl_, &
            e1,o1
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='read', access='direct',recl=recl)
       read(un,rec=ikr) nstval_,nstcon_,nkpt_,ngq_,vql_,vkl_, &
            e1,o1
    end if
    close(un)
    ! check consistency
    if ((nstval_.ne.nstval).or.(nstcon_.ne.nstcon).or.(nkpt_.ne.nkpt).or. &
         ( r3dist(vql_,vql(1,iq)) > epslat ).or. &
         ( r3dist(vkl_,vkl(1,ik)) > epslat )) then
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &matrix elements (current/file): '
       write(unitout,'(a,2i6)') 'nstval', nstval, nstval_
       write(unitout,'(a,2i6)') 'nstcon', nstcon, nstcon_
       write(unitout,'(a,2i6)') 'nkpt', nkpt, nkpt_
       write(unitout,'(a,3f12.6,a,3f12.6)') 'vql', vql(:,iq), ',', vql_
       write(unitout,'(a,3f12.6,a,3f12.6)') 'vkl', vkl(:,ik), ',', vkl_
       write(unitout,'(a)') ' file: ',trim(filnam)
       call terminate
    end if

  end subroutine getdevalsv2

end module m_getdevalsv2
