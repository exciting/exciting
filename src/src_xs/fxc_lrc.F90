


! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_fxc_lrc
  implicit none

contains

!BOP
! !ROUTINE: fxc_lrc
! !INTERFACE:


subroutine fxc_lrc(msiz, sw, alpha, fxc)
! !USES:
    use modmain
    use modxs
! !INPUT/OUTPUT PARAMETERS:
!   msiz  : matrix size of local field effects (in,integer)
!   sw    : true for inclusion of local field effects (in,logical)
!   alpha : real constant (in,real)
!   fxc   : xc-kernel Fourier coefficients (out,complex(:,:))
! !DESCRIPTION:
!   Static long range xc-kernel; S. Botti, PRB 70, 045301 (2004).
!   Calculates the symmetrized xc-kernel for the static long range model.
!   According to the switch {\tt sw} either $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{\alpha}{4\pi} 
!         \delta({\bf G},{\bf G'}), $$
!   if the switch is true, or $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{\alpha}{4\pi} 
!         \delta({\bf G},{\bf G'})\delta({\bf G},{\bf 0}), $$
!   otherwise.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: msiz
    ! true if all G-components of fxc are to be considered
    logical, intent(in) :: sw
    real(8), intent(in) :: alpha
    complex(8), intent(out) :: fxc(:, :)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_lrc'
    real(8) :: t1
    integer :: sh(2), ig
    sh=shape(fxc)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout, '(a, 2i9, a, i9, a)') 'Error('//trim(thisnam)//'): size of &
	    &fxc is to small (required)', sh, '(', msiz, ')'
       call terminate
    end if
    fxc(:, :)=(0.d0, 0.d0)
    t1=-alpha/fourpi
    if (.not.sw) then
       fxc(1, 1)=t1
    else
       do ig=1, msiz
	  fxc(ig, ig)=t1
       end do
    end if
  end subroutine fxc_lrc
end module m_fxc_lrc
!EOC
