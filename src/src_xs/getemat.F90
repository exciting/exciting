
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getemat
  implicit none
contains

  subroutine getemat(iq,ik,tarec,n1,n2,n3,n4,ngp,filnam,x1,x2)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,n1,n2,n3,n4,ngp
    logical :: tarec
    character(*) :: filnam
    complex(8), intent(out) :: x1(:,:,:)
    complex(8), optional, intent(out) :: x2(:,:,:)
    ! local variables
    character(*), parameter :: thisnam = 'getemat'
    integer :: recl,un,ikr,n1_,n2_,n3_,n4_,nstsv_,ngq_,err
    real(8) :: vql_(3),vkl_(3)
    logical :: existent, lerr
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
    ! check block sizes
    lerr=(n1.lt.1).or.(n1.gt.nstsv).or.(n2.lt.1).or.(n2.gt.nstsv)
    if (present(x2)) lerr=lerr.or.(n3.lt.1).or.(n3.gt.nstsv).or.(n4.lt.1).or. &
         (n4.gt.nstsv)
    if (lerr) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): inconsistent limits for bands:")') &
            thisnam
       if (present(x2)) then
          write(unitout,'(" requested band block sizes : ",4i6)') n1,n2,n3,n4
       else
          write(unitout,'(" requested band block sizes : ",4i6)') n1,n2
       end if
       write(unitout,'(" maximum value              : ",i6)') nstsv
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    lerr=(size(x1,1).ne.n1).or.(size(x1,2).ne.n2).or.(ngp.gt.ngq(iq))
    if (present(x2)) lerr=lerr.or.(size(x2,1).ne.n3).or.(size(x2,2).ne.n4)
    if (lerr) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): output array does not match for &
            &bands:")') thisnam
       write(unitout,'(" requested number of G+q vectors : ",i6)') ngp
       write(unitout,'(" current number of G+q vectors   : ",i6)') ngq(iq)
       if (present(x2)) then
          write(unitout,'(" block sizes : ",4i6)') n1,n2,n3,n4
          write(unitout,'(" array sizes : ",4i6)') size(x1,1),size(x1,2), &
               size(x2,1),size(x2,2)
       else
          write(unitout,'(" block sizes : ",2i6)') n1,n2
          write(unitout,'(" array sizes : ",2i6)') size(x1,1),size(x1,2)
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
    if (present(x2)) then
       ! I/O record length
       inquire(iolength=recl) vql_,vkl_,nstsv_,ngq_,n1_,n2_,n3_,n4_, &
            x1,x2
       open(unit=un,file=trim(filnam),form='unformatted',action='read', &
            access='direct',recl=recl)
       read(un,rec=ikr) vql_,vkl_,nstsv_,ngq_,n1_,n2_,n3_,n4_,x1,x2
    else
       ! I/O record length
       inquire(iolength=recl) vql_,vkl_,nstsv_,ngq_,n1_,n2_,x1
       open(unit=un,file=trim(filnam),form='unformatted',action='read', &
            access='direct',recl=recl)
       read(un,rec=ikr) vql_,vkl_,nstsv_,ngq_,n1_,n2_,x1
    end if
    close(un)





!!$    ! check consistency
!!$    if ((n1_.ne.n1).or.(n2_.ne.n2).or. &
!!$         (r3dist(vql_,vql(1,iq)).gt.epslat).or. &
!!$         (r3dist(vkl_,vkl(1,ik)).gt.epslat)) then
!!$       write(unitout,*)
!!$       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
!!$            &matrix elements (current/file): '
!!$       write(unitout,'(a,2i6)') ' n1:', n1, n1_
!!$       write(unitout,'(a,2i6)') ' n2:', n2, n2_
!!$       if (present(x2)) then
!!$          write(unitout,'(a,2i6)') ' n3:', n3, n3_
!!$          write(unitout,'(a,2i6)') ' n4:', n4, n4_
!!$       end if
!!$       write(unitout,'(a,3f12.6,a,3f12.6)') ' vql :', vql(:,iq), ',', vql_
!!$       write(unitout,'(a,3f12.6,a,3f12.6)') ' vkl :', vkl(:,ik), ',', vkl_
!!$       write(unitout,'(a)') ' file: ',trim(filnam)
!!$       write(unitout,*)
!!$       call flushifc(unitout)
!!$       call terminate
!!$    end if



  end subroutine getemat

end module m_getemat
