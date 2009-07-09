


! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putx0
  implicit none
contains


subroutine putx0(tp0, iq, iw, filnam, filxt, ch0, ch0wg, ch0hd)
    use modmain
    use modxs
    use m_getunit
    implicit none
    ! arguments
    logical, intent(in) :: tp0
    integer, intent(in) :: iq, iw
    character(*), intent(in) :: filnam, filxt
    complex(8), intent(in) :: ch0(:, :)
    complex(8), intent(in), optional :: ch0wg(:, :, :), ch0hd(:, :)
    ! local variables
    character(*), parameter :: thisnam = 'putx0'
    integer :: un, recl
    ! q=0 but head or wings missing
    if (tp0.and.((.not.present(ch0wg)).or.(.not.present(ch0wg))) ) then
       write(*, *) 'Error('//trim(thisnam)//'): q=0 but head or wings missing'
       call terminate
    end if
    call getunit(un)
    if (tp0) then
       ! I/O record length
       inquire(iolength=recl) ngq(iq), vql(:, iq), ch0, ch0wg, ch0hd
       open(unit = un, file = trim(filnam)//trim(filxt), form = 'unformatted', &
	    action = 'write', access = 'direct', recl = recl)
       write(un, rec=iw) ngq(iq), vql(:, iq), ch0, ch0wg, ch0hd
    else
       ! I/O record length
       inquire(iolength=recl) ngq(iq), vql(:, iq), ch0
       open(unit = un, file = trim(filnam)//trim(filxt), form = 'unformatted', &
	    action = 'write', access = 'direct', recl = recl)
       write(un, rec=iw) ngq(iq), vql(:, iq), ch0
    end if
    close(un)
  end subroutine putx0

end module m_putx0
