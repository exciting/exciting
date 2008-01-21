
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getemat2
  implicit none
contains

  subroutine getemat2(iq,ik,tarec,filnam,x1,x2)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    logical :: tarec
    character(*) :: filnam
    complex(8), intent(out) :: x1(:,:,:)
    complex(8), optional, intent(out) :: x2(:,:,:)
    ! local variables
    character(*), parameter :: thisnam = 'getemat2'
    integer :: recl,un,ikr,nst1_,nst2_,nst3_,nst4_,nkpt_,ngq_
    real(8) :: vql_(3),vkl_(3)
    logical :: existent
    ! functions
    real(8) :: r3dist
    external :: r3dist
    ! check if file exists
    inquire(file=trim(filnam),exist=existent)
    if (.not.existent) then
       write(unitout,'(a)') 'Error('//thisnam//'): file does not exist: '// &
            trim(filnam)
       call terminate
    end if
    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)
    call getunit(un)
    if (present(x2)) then
       ! I/O record length
       inquire(iolength=recl) nst1_,nst2_,nst3_,nst4_,nkpt_,ngq_,vql_,vkl_, &
            x1,x2
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='read', access='direct',recl=recl)
       read(un,rec=ikr) nst1_,nst2_,nst3_,nst4_,nkpt_,ngq_,vql_,vkl_,x1,x2
    else
       ! I/O record length
       inquire(iolength=recl) nst1_,nst2_,nkpt_,ngq_,vql_,vkl_,x1
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='read', access='direct',recl=recl)
       read(un,rec=ikr) nst1_,nst2_,nkpt_,ngq_,vql_,vkl_,x1
    end if
    close(un)
    ! check consistency
    if ((nst1_.ne.nst1).or.(nst2_.ne.nst2).or. &
         (nkpt_.ne.nkpt).or. &
         (r3dist(vql_,vql(1,iq)).gt.epslat).or. &
         (r3dist(vkl_,vkl(1,ik)).gt.epslat)) then
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &matrix elements (current/file): '
       write(unitout,'(a,2i6)') ' nst1:', nst1, nst1_
       write(unitout,'(a,2i6)') ' nst2:', nst2, nst2_
       if (present(x2)) then
          write(unitout,'(a,2i6)') ' nst3:', nst3, nst3_
          write(unitout,'(a,2i6)') ' nst4:', nst4, nst4_
       end if
       write(unitout,'(a,2i6)') ' nkpt', nkpt, nkpt_
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vql :', vql(:,iq), ',', vql_
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vkl :', vkl(:,ik), ',', vkl_
       write(unitout,'(a)') ' file: ',trim(filnam)
       call terminate
    end if
  end subroutine getemat2

end module m_getemat2
