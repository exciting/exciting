

! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_puttetcw
  implicit none
contains


subroutine puttetcw(iq, ik, i1, i2, n1, n2, filnam, cw, cwa, cwsurf)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq, ik, i1, i2, n1, n2
    character(*), intent(in) :: filnam
    real(8), intent(in) :: cw(:), cwa(:), cwsurf(:)
    ! local variables
    integer :: un, recl, irec, iqt, err
    ! check input parameters
    err=0
    if ((n1.lt.1).or.(n1.gt.nstsv)) then
       write(unitout, *)
       write(unitout, '("Error(puttetcw): n1 < 1 or n1 > nstsv")')
       write(unitout, '(" n1	:", i6)') n1
       write(unitout, '(" nstsv :", i6)') nstsv
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    if ((n2.lt.1).or.(n2.gt.nstsv)) then
       write(unitout, *)
       write(unitout, '("Error(puttetcw): n2 < 1 or n2 > nstsv")')
       write(unitout, '(" n1	:", i6)') n2
       write(unitout, '(" nstsv :", i6)') nstsv
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.ne.0) call terminate
    ! q-point
    iqt=iq
    ! record position
    irec=(ik-1)*n1*n2 + (i1-1)*n2 + i2
    ! I/O record length
    inquire(iolength=recl) vql(:, iq), vkl(:, ik), nstsv, n1, n2, cw, cwa, cwsurf
    call getunit(un)
    open(unit = un, file = trim(filnam), form = 'unformatted', &
	 action = 'write', access = 'direct', recl = recl)
    write(un, rec=irec) vql(:, iq), vkl(:, ik), nstsv, n1, n2, cw, cwa, cwsurf
    close(un)
  end subroutine puttetcw

end module m_puttetcw
