
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getemat
  implicit none
contains

  subroutine getemat(iq,ik,tarec,filnam,ngp,l1,h1,l2,h2,x12,l3,h3,l4,h4,x34)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,ngp
    logical :: tarec
    character(*) :: filnam
    integer, intent(in) :: l1,h1,l2,h2
    complex(8), intent(out) :: x12(:,:,:)
    integer, optional, intent(in) :: l3,h3,l4,h4
    complex(8), optional, intent(out) :: x34(:,:,:)
    ! local variables
    character(*), parameter :: thisnam = 'getemat'
    integer :: recl,un,ikr,n1,n2,n3,n4,nstsv_,ngq_,err
    integer :: l1_,h1_,l2_,h2_,l3_,h3_,l4_,h4_,n1_,n2_,n3_,n4_
    real(8) :: vql_(3),vkl_(3)
    logical :: existent, lerr
    complex(8), allocatable :: x12t(:,:,:),x34t(:,:,:)
    ! functions
    real(8) :: r3dist
    external :: r3dist
    ! check if all optional variables are present if any is present
    if ((present(l3).or.present(h3).or.present(l4).or.present(h4).or. &
         present(x34)).and.(.not.( &
         present(l3).and.present(h3).and.present(l4).and.present(h4).and. &
         present(x34)))) then
       write(*,*)
       write(*,'("Error(getemat): optional parameters not complete - check &
            &calling routines")')
       write(*,*)
       call terminate
    end if
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
    ! check limits for states
    lerr=(l1.lt.1).or.(l1.gt.nstsv).or.(h1.lt.1).or.(h1.gt.nstsv).or. &
         (l2.lt.1).or.(l2.gt.nstsv).or.(h2.lt.1).or.(h2.gt.nstsv).or. &
         (l1.ge.h1).or.(l2.ge.h2)
    if (present(x34)) &
         lerr=lerr.or.(l3.lt.1).or.(l3.gt.nstsv).or.(h3.lt.1).or.(h3.gt.nstsv) &
         .or.(l4.lt.1).or.(l4.gt.nstsv).or.(h4.lt.1).or.(h4.gt.nstsv).or. &
         (l3.ge.h3).or.(l4.ge.h4)
    err=0
    if (lerr) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): inconsistent requested limits for &
            &states:")') thisnam
       if (present(x34)) then
          write(unitout,'(" requested state limits (lo,hi): ",4(2i6,2x))') &
               l1,h1,l2,h2,l3,h3,l4,h4
       else
          write(unitout,'(" requested state limits (lo,hi): ",2(2i6,2x))') &
               l1,h1,l2,h2
       end if
       write(unitout,'(" maximum value                 : ",i6)') nstsv
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    n1=h1-l1+1
    n2=h2-l2+1
    if (present(x34)) then
       n3=h3-l3+1
       n4=h4-l4+1
    end if
    ! check block sizes against array
    lerr=(size(x12,1).ne.n1).or.(size(x12,2).ne.n2).or.(ngp.gt.ngq(iq)) &
         .or.(size(x12,3).ne.ngp)
    if (present(x34)) lerr=lerr.or.(size(x34,1).ne.n3).or.(size(x34,2).ne.n4) &
         .or.(size(x34,3).ne.ngp)
    if (lerr) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): output array does not match for &
            &states:")') thisnam
       write(unitout,'(" requested number of G+q vectors : ",i6)') ngp
       write(unitout,'(" current number of G+q vectors   : ",i6)') ngq(iq)
       if (present(x34)) then
          write(unitout,'(" array sizes for G+q vectors     : ",2i6)') &
               size(x12,3),size(x34,3)
          write(unitout,'(" block sizes : ",4i6)') n1,n2,n3,n4
          write(unitout,'(" array sizes : ",4i6)') size(x12,1),size(x12,2), &
               size(x34,1),size(x34,2)
       else
          write(unitout,'(" array size for G+q vectors      : ",i6)') &
               size(x12,3)
          write(unitout,'(" block sizes : ",2i6)') n1,n2
          write(unitout,'(" array sizes : ",2i6)') size(x12,1),size(x12,2)
       end if
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
    call getunit(un)
    if (present(x34)) then
       ! I/O record length
       inquire(iolength=recl) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_,l3_,h3_, &
            l4_,h4_
       open(unit=un,file=trim(filnam),form='unformatted',action='read', &
            access='direct',recl=recl)
       read(un,rec=1) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_,l3_,h3_,l4_,h4_
       close(un)
    else
       ! I/O record length
       inquire(iolength=recl) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_
       open(unit=un,file=trim(filnam),form='unformatted',action='read', &
            access='direct',recl=recl)
       read(un,rec=1) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_
       close(un)
    end if
    err=0
    ! check block sizes
    lerr=(l1.lt.l1_).or.(h1.gt.h1_).or.(l2.lt.l2_).or.(h2.gt.h2_)
    if (present(x34)) lerr=lerr.or.(l3.lt.l3_).or.(h3.gt.h3_).or. &
         (l4.lt.l4_).or.(h4.gt.h4_)
    if (lerr) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): limits for states out of range in &
            &file:")') thisnam
       if (present(x34)) then
          write(unitout,'(" requested state limits (lo,hi): ",4(2i6,2x))') &
               l1,h1,l2,h2,l3,h3,l4,h4
          write(unitout,'(" state limits from file (lo,hi): ",4(2i6,2x))') &
               l1_,h1_,l2_,h2_,l3_,h3_,l4_,h4_
       else
          write(unitout,'(" requested state limits (lo,hi): ",2(2i6,2x))') &
               l1,h1,l2,h2
          write(unitout,'(" state limits from file (lo,hi): ",2(2i6,2x))') &
               l1_,h1_,l2_,h2_
          write(unitout,'(" file                          : ",a)') trim(filnam)
       end if
       write(unitout,*)
       call flushifc(unitout)
       err=err+1 
    end if
    if (err.gt.0) call terminate
    !------------------!
    !     get data     !
    !------------------!
    n1_=h1_-l1_+1; n2_=h2_-l2_+1
    if (present(x34)) then
       n3_=h3_-l3_+1; n4_=h4_-l4_+1
    end if
    allocate(x12t(n1_,n2_,ngq_))
    if (present(x34)) allocate(x34t(n3_,n4_,ngq_))
    call getunit(un)
    if (present(x34)) then
       ! I/O record length
       inquire(iolength=recl) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_,l3_,h3_, &
            l4_,h4_,x12t,x34t
       open(unit=un,file=trim(filnam),form='unformatted',action='read', &
            access='direct',recl=recl)
       read(un,rec=ikr) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_,l3_,h3_,l4_,h4_,&
            x12t,x34t
    else
       ! I/O record length
       inquire(iolength=recl) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_,x12t
       open(unit=un,file=trim(filnam),form='unformatted',action='read', &
            access='direct',recl=recl)
       read(un,rec=ikr) vql_,vkl_,nstsv_,ngq_,l1_,h1_,l2_,h2_,x12t
    end if
    close(un)
    ! check q-point and k-point
    if ((r3dist(vql_,vql(1,iq)).gt.epslat).or. &
         (r3dist(vkl_,vkl(1,ik)).gt.epslat)) then
       write(unitout,*)
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &matrix elements (current/file): '
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vql :', vql(:,iq), ',', vql_
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vkl :', vkl(:,ik), ',', vkl_
       write(unitout,'(a)') ' file: ',trim(filnam)
       write(unitout,*)
       call flushifc(unitout)
       call terminate
    end if
    ! retrieve data within cutoff
    x12(:,:,:)=x12t(l1-l1_+1:h1-l1_+1,l2-l2_+1:h2-l2_+1,:ngp)
    deallocate(x12t)
    if (present(x34)) x34(:,:,:)=x34t(l3-l3_+1:h3-l3_+1,l4-l4_+1:h4-l4_+1,:ngp)
    if (present(x34)) deallocate(x34t)
  end subroutine getemat

end module m_getemat
