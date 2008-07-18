
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gettetcw
  implicit none
contains
  
  subroutine gettetcw(iq,ik,i1,i2,n1,n2,nw,fnam,cw,cwa,cwsurf)
    use modmain
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,i1,i2,n1,n2,nw
    character(*), intent(in) :: fnam
    real(8), intent(out) :: cw(nw),cwa(nw),cwsurf(nw)
    ! local variables
    integer :: un,irec,recl,nstsv_,n1_,n2_,err
    real(8) :: vql_(3),vkl_(3),vqlt(3),vklt(3)
    real(8), external :: r3dist
    err=0
    if ((i1.lt.1).or.(i1.gt.nstsv).or.(i2.lt.1).or.(i2.gt.nstsv)) then
       write(unitout,*)
       write(unitout,'("Error(gettetcw): inconsistent band combination:")')
       write(unitout,'(" bands         : ",2i6)') i1,i2
       write(unitout,'(" maximum value : ",i6)') nstfv
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
    inquire(iolength=recl) vql_,vkl_,nstsv_,n1_,n2_
    call getunit(un)
    open(un,file=trim(fnam),action='read',form='unformatted',status='old', &
         access='direct',recl=recl)
    read(un,rec=1) vql_,vkl_,nstsv_,n1_,n2_
    close(un)
    err=0
    ! check number of bands
    if (nstsv.gt.nstsv_) then
       write(unitout,*)
       write(unitout,'("Error(gettetcw): invalid nstsv for k-point ",I8)') ik
       write(unitout,'(" q-point    : ",I8)') iq
       write(unitout,'(" current    : ",I8)') nstsv
       write(unitout,'(" FILE       : ",I8)') nstsv_
       write(unitout,'(" filename   : ",a )') trim(fnam)
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if ((n1.ne.n1_).or.(n2.gt.n2_)) then
       write(unitout,*)
       write(unitout,'("Error(gettetcw): invalid band ranges for k-point ",&
            &I8)') ik
       write(unitout,'(" q-point    : ",I8)') iq
       write(unitout,'(" current    : ",2I8)') n1,n2
       write(unitout,'(" FILE       : ",2I8)') n1_,n2_
       write(unitout,'(" filename   : ",a )') trim(fnam)
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! record position with proper n1 and n2 values
    irec=(ik-1)*n1_*n2_ + (i1-1)*n2_ + i2
    ! read from file
    call getunit(un)
    inquire(iolength=recl) vql_,vkl_,nstsv_,n1_,n2_,cw,cwa,cwsurf
    open(un,file=trim(fnam),form='unformatted',action='read',status='old', &
         access='direct',recl=recl)
    read(un,rec=irec) vql_,vkl_,nstsv_,n1_,n2_,cw,cwa,cwsurf
    close(un)
    ! check q-point and k-point
    if (iq.eq.0) then
       ! Gamma Q-point
       vklt(:)=vkl0(:,ik)
       vqlt(:)=0.d0
    else
       vklt(:)=vkl(:,ik)
       vqlt(:)=vql(:,iq)
    end if
    if ((r3dist(vkl_,vklt).gt.epslat).or.(r3dist(vql_,vqlt).gt.epslat)) then
       write(unitout,*)
       write(unitout,'(a)') 'Error(gettetcw): differring parameters for &
            &tetrahedron convolution weights (current/file): '
       write(unitout,'(a,i6)') ' q-point index  :', iq
       write(unitout,'(a,i6)') ' k-point index  :', ik
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vql            :', vqlt,',', vql_
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vkl            :', vklt,',', vkl_
       write(unitout,'(a)')    ' file           : '//trim(fnam)
       write(unitout,*)
       call flushifc(unitout)
       call terminate
    end if
  end subroutine gettetcw

end module m_gettetcw
