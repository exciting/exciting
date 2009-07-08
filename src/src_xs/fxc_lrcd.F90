

! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_fxc_lrcd
  implicit none

contains

!BOP
! !ROUTINE: fxc_lrcd
! !INTERFACE:


subroutine fxc_lrcd(msiz, sw, alpha, beta, w, fxc)
! !USES:
    use modmain
    use modxs
! !INPUT/OUTPUT PARAMETERS:
!   msiz  : matrix size of local field effects (in,integer)
!   sw    : true for inclusion of local field effects (in,logical)
!   alpha : real constant (in,real)
!   w     : frequency grid (in,complex(:))
!   fxc   : xc-kernel Fourier coefficients (out,complex(:,:))
! !DESCRIPTION:
!   Dynamical long range xc-kernel; S. Botti, PRB 72, 125203 (2005).
!   Calculates the symmetrized xc-kernel for the static long range model.
!   According to the switch {\tt sw} either $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{1}{4\pi} (\alpha+\beta\omega^2)
!         \delta({\bf G},{\bf G'}), $$
!   if the switch is true, or $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{1}{4\pi} (\alpha+\beta\omega^2)
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
    real(8), intent(in) :: alpha, beta
    complex(8), intent(in) :: w
    complex(8), intent(out) :: fxc(:, :)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_lrcd'
    complex(8) :: zt1
    integer :: sh(2), ig
    sh=shape(fxc)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout, '(a, 2i9, a, i9, a)') 'Error('//trim(thisnam)//'): size of &
	    &fxc is to small (required)', sh, '(', msiz, ')'
       call terminate
    end if    
    fxc(:, :)=(0.d0, 0.d0)
    zt1=-(alpha+beta*w**2)/fourpi
    if (.not.sw) then
       fxc(1, 1)=zt1
    else
       do ig=1, msiz
	  fxc(ig, ig)=zt1
       end do
    end if
  end subroutine fxc_lrcd
!EOC

end module m_fxc_lrcd
