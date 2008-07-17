
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getpmat
  implicit none
contains

  subroutine getpmat(ik,vklt,isti,istf,tarec,filnam,pm)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: ik,isti,istf
    real(8), intent(in) :: vklt(:,:)
    logical :: tarec
    character(*) :: filnam
    complex(8), intent(out) :: pm(:,:,:)
    ! local variables
    character(*), parameter :: thisnam='getpmat'
    integer :: recl,un,ikr,nstsv_,err
    real(8) :: vkl_(3)
    logical :: existent
    complex(8), allocatable :: pmt(:,:,:)
    ! functions
    real(8), external :: r3dist
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
    err=0
    ! check band range
    if ((isti.lt.1).or.(istf.gt.nstsv).or.(istf.le.isti)) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): inconsistent limits for bands:")') thisnam
       write(unitout,'(" band limits   : ",2i6)') isti,istf
       write(unitout,'(" maximum value : ",i6)') nstsv
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if ((size(pm,2).ne.(istf-isti+1)).or.(size(pm,3).ne.(istf-isti+1))) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): output array does not match for &
            &bands:")') thisnam
       write(unitout,'(" band limits               : ",2i6)') isti,istf
       write(unitout,'(" requested number of bands : ",i6)') istf-isti+1
       write(unitout,'(" array sizes               : ",2i6)') size(pm,2), &
            size(pm,3)
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
    inquire(iolength=recl) vkl_,nstsv_
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted',action='read', &
         access='direct',recl=recl)
    read(un,rec=1) vkl_,nstsv_
    close(un)
    err=0
    ! check number of bands
    if (nstsv.gt.nstsv_) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): invalid nstsv for k-point ",I8)') &
            thisnam,ik
       write(unitout,'(" current    : ",I8)') nstsv
       write(unitout,'(" FILE       : ",I8)') nstsv_
       write(unitout,'(" filename   : ",a )') trim(filnam)
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! allocate local arrays
    allocate(pmt(3,nstsv_,nstsv_))
    ! I/O record length
    inquire(iolength=recl) vkl_,nstsv_,pmt
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted',action='read', &
         access='direct',recl=recl)
    ! read from file
    read(un,rec=ikr) vkl_,nstsv_,pmt
    close(un)
    ! check k-point
    if (r3dist(vkl_,vklt(1,ik)).gt.epslat) then
       write(unitout,*)
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &matrix elements (current/file): '
       write(unitout,'(a,i6)') ' k-point index  :', ik
       write(unitout,'(a,i6)') ' record position:', ikr
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vkl            :', vklt(:,ik), &
            ',', vkl_
       write(unitout,'(a)')    ' file           : ',trim(filnam)
       write(unitout,*)
       call flushifc(unitout)
       call terminate
    end if
    ! retrieve data within cutoffs
    pm(:,:,:)=pmt(:,isti:istf,isti:istf)
    deallocate(pmt)
  end subroutine getpmat

end module m_getpmat
