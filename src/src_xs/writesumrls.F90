

! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writesumrls
  implicit none
contains


subroutine writesumrls(iq, s, fn)
    use modmain
    use modxs
    use m_getunit
    use m_writevars
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: s(3)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam = 'writesumrls'
    call getunit(unit1)
    open(unit1, file=trim(fn), action='write')
    ! zeroth frequency moment sumrule
    write(unit1, '(a, g18.10, a, g18.10, a)') 'zeroth frequency moment sumrule &
	 &(num. val. el.):', s(1), '(', chgval/2.d0, ')'
    ! first frequency moment sumrule
    write(unit1, '(a, g18.10, a, g18.10, a)') 'first frequency moment sumrule  &
	 &(num. val. el.):', s(2), '(', chgval/2.d0, ')'
    ! one over frequency sumrule
    write(unit1, '(a, g18.10, a, g18.10, a)') 'pi half sumrule		       &
	 &(target)	 :', s(3), '(', pi/2.d0, ')'
    ! write parameters as header to file
    call writevars(unit1, iq, iq)
    ! close file
    close(unit1)
  end subroutine writesumrls

end module m_writesumrls
