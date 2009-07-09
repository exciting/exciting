


! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getpmat
  implicit none
contains


subroutine getpmat(ik, vklt, i1, f1, i2, f2, tarec, filnam, pm)
    use modmain
use modinput
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: ik, i1, f1, i2, f2
    real(8), intent(in) :: vklt(:, :)
    logical :: tarec
    character(*) :: filnam
    complex(8), intent(out) :: pm(:, :, :)
    ! local variables
    character(*), parameter :: thisnam='getpmat'
    integer :: recl, un, ikr, nstsv_, err
    real(8) :: vkl_(3)
    logical :: existent
    complex(8), allocatable :: pmt(:, :, :)
    ! functions
    real(8), external :: r3dist
    ! check if file exists
    inquire(file=trim(filnam), exist=existent)
    if (.not.existent) then
       write(unitout, '(a)') 'Error('//thisnam//'): file does not exist: '// &
	    trim(filnam)
       call terminate
    end if
    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(procs, nkpt, ik, ikr)
    err=0
    ! check band range
    if ((i1.lt.1).or.(i1.gt.nstsv).or.(f1.lt.1).or.(f1.gt.nstsv).or. &
	 (i2.lt.1).or.(i2.gt.nstsv).or.(f2.lt.1).or.(f2.gt.nstsv).or. &
	 (i1.ge.f1).or.(i2.ge.f2)) then
       write(unitout, *)
       write(unitout, '("Error(", a, "): inconsistent limits for states:")') &
	    thisnam
       write(unitout, '(" limits (lo/hi) : ", 2(2i6, 2x))') i1, f1, i2, f2
       write(unitout, '(" maximum value  : ", i6)') nstsv
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    if ((size(pm, 2).ne.(f1-i1+1)).or.(size(pm, 3).ne.(f2-i2+1))) then
       write(unitout, *)
       write(unitout, '("Error(", a, "): output array does not match for &
	    &states:")') thisnam
       write(unitout, '(" limits		     : ", 2(2i6, 2x))') i1, f1, i2, f2
       write(unitout, '(" requested number of states : ", 2i6)') f1-i1+1, f2-i2+1
       write(unitout, '(" array sizes		     : ", 2i6)') size(pm, 2), &
	    size(pm, 3)
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
    inquire(iolength=recl) vkl_, nstsv_
    call getunit(un)
    open(unit = un, file = trim(filnam), form = 'unformatted', action = 'read', &
	 access = 'direct', recl = recl)
    read(un, rec=1) vkl_, nstsv_
    close(un)
    err=0
    ! check if all states can be read from file
    if ((f1.gt.nstsv_).or.(f2.gt.nstsv_)) then
       write(unitout, *)
       write(unitout, '("Error(", a, "): requested states out of range for &
	    &k-point ", I8)') thisnam, ik
       write(unitout, '(" limits	  : ", 2(2I6, 2x))') i1, f1, i2, f2
       write(unitout, '(" cutoff from file: ", I8)') nstsv_
       write(unitout, '(" filename	  : ", a )') trim(filnam)
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! allocate local arrays
    allocate(pmt(3, nstsv_, nstsv_))
    ! I/O record length
    inquire(iolength=recl) vkl_, nstsv_, pmt
    call getunit(un)
    open(unit = un, file = trim(filnam), form = 'unformatted', action = 'read', &
	 access = 'direct', recl = recl)
    ! read from file
    read(un, rec=ikr) vkl_, nstsv_, pmt
    close(un)
    ! check k-point
    if (r3dist(vkl_, vklt(1, ik)).gt.input%structure%epslat) then
       write(unitout, *)
       write(unitout, '(a)') 'Error('//thisnam//'): differring parameters for &
	    &matrix elements (current/file): '
       write(unitout, '(a, i6)') ' k-point index  :', ik
       write(unitout, '(a, i6)') ' record position:', ikr
       write(unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :', vklt(:, ik), &
	    ', ', vkl_
       write(unitout, '(" filename	  : ", a )') trim(filnam)
       write(unitout, *)
       call flushifc(unitout)
       call terminate
    end if
    ! retrieve data within cutoffs
    pm(:, :, :)=pmt(:, i1:f1, i2:f2)
    deallocate(pmt)
  end subroutine getpmat

end module m_getpmat
