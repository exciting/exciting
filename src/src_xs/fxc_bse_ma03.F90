



! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_fxc_bse_ma03
  implicit none

contains

!BOP
! !ROUTINE: fxc_bse_ma03
! !INTERFACE:


subroutine fxc_bse_ma03(msiz, oct, sw, iw, fxc)
! !USES:
use modinput
    use mod_constants, only:zzero
    use modmpi, only:
    use modxs, only:unitout, bzsampl
    use invert
    use m_xsgauntgen
    use m_findgntn0
    use m_writegqpts
    use m_genfilname
    use m_getunit
! !INPUT/OUTPUT PARAMETERS:
!   msiz  : matrix size of local field effects (in,integer)
!   sw    : true for inclusion of local field effects (in,logical)
!   alpha : real constant (in,real)
!   fxc   : xc-kernel Fourier coefficients (out,complex(:,:))
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003).
!   Interface function.
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: msiz, oct
    ! true if all G-components of fxc are to be considered
    logical, intent(in) :: sw
    integer, intent(in) :: iw
    complex(8), intent(out) :: fxc(:, :)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_bse_ma03'
    character(256) :: filnam
    complex(8), allocatable :: fxch(:, :), fxcw1(:, :), fxcw2(:, :)
    complex(8) :: zt1
    integer :: n, n2, un, recl, n_
    n=size(fxc, 1)
    n2=size(fxc, 2)
    if ((n.lt.msiz).or.(n.ne.n2)) then
       write(unitout, '(a, 2i9, a, i9, a)') 'Error('//trim(thisnam)//'): size of &
	    &fxc is inconsistent (required)', n, n2, '(', msiz, ')'
       call terminate
    end if
    allocate(fxch(-3:-1, -3:-1), fxcw1(-3:-1, n), fxcw2(n, -3:-1))
    ! filename for BSE-xc-kernel
    call getunit(un)
    ! filename for xc-kernel
    call genfilname(basename = 'FXC_BSE', asc = .false., bzsampl = bzsampl, &
	 acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresfxc, tord=input%xs%tddft%tordfxc, iqmt = 1, filnam&
    &= filnam)
	  ! get LFE size
    inquire(iolength=recl) n_
    open(un, file = trim(filnam), form = 'unformatted', action = 'read', &
	 status = 'old', access = 'direct', recl = recl)
    read(un, rec=1) n_
    close(un)
    ! check if data from file can be stored in local array
    if (n.lt.n_) then
       write(unitout, *)
       write(unitout, '("Error(", a, "): LFE size of file too large")') trim(thisnam)
       write(unitout, '(" LFE size (file)    : ", i8)') n_
       write(unitout, '(" LFE size (current) : ", i8)') n
       write(unitout, *)
       call terminate
    end if
	  ! get data
	  fxc(:, :)=zzero
    inquire(iolength=recl) n_, fxch, fxcw1(:, :n_), fxcw2(:n_, :), fxc(:n_, :n_)
    open(un, file = trim(filnam), form = 'unformatted', action = 'read', &
	 status = 'old', access = 'direct', recl = recl)
    read(un, rec=iw) n_, fxch, fxcw1(:, :n_), fxcw2(:n_, :), fxc(:n_, :n_)
    ! assign head
    fxc(1, 1)=fxch(-oct, -oct)
    ! assign wings
    if (msiz.gt.1) then
       fxc(1, 2:n_)=fxcw1(-oct, 2:n_)
       fxc(2:n_, 1)=fxcw2(2:n_, -oct)
    end if
    ! no LFE at all
    if (.not.sw) then
       zt1=fxc(1, 1)
       fxc(:, :)=zzero
       fxc(1, 1)=zt1
    end if
    close(un)
    deallocate(fxch, fxcw1, fxcw2)
  end subroutine fxc_bse_ma03
!EOC

end module m_fxc_bse_ma03
