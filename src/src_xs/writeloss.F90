


! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writeloss
  implicit none
contains


subroutine writeloss(iq, w, loss, fn)

    use mod_lattice
    use mod_constants
    use mod_charge_and_moment
    use modxs
    use m_getunit
    use m_writevars
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: w(:)
    real(8), intent(in) :: loss(:)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam='writeloss'
    integer :: n1(1), n, iw, igmt
    if (any(shape(w).ne.shape(loss))) then
       write(unitout, '(a)') 'Error('//thisnam//'): input arrays have &
	    &diffenrent shape'
       call terminate
    end if
    n1=shape(w); n=n1(1)
    call getunit(unit1)
    open(unit1, file=trim(fn), action='write')
    ! include dynamical structure factor
    ! Dynamical structure factor; expression taken from Weissker, PRL 2006
    ! Units of dynamical structure factor are Hartree^-1
    igmt=ivgigq(ivgmt(1, iq), ivgmt(2, iq), ivgmt(3, iq), iq)
    write(unit1, '(3g18.10)') (w(iw) * escale, loss(iw), loss(iw)* &
	 (gqc(igmt, iq) ** 2/(4.d0 * pi ** 2 * chgval/omega)), iw = 1, n)
    ! write relevant parameters to file
    call writevars(unit1, iq, iq)
    close(unit1)
  end subroutine writeloss

end module m_writeloss
