
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gensigma
  implicit none
contains

  subroutine gensigma(w,eps,oc,sigma)
    use modmain
    use modxs
    implicit none
    ! arguments
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: eps(:)
    integer, intent(in) :: oc(3)
    complex(8), intent(out) :: sigma(:)
    ! local variables
    character(*), parameter :: thisnam = 'gensigma'
    real(8) :: delt

    if (any(shape(eps).ne.shape(sigma))) then
       write(unitout,'(a)') 'Error('//thisnam//'): input and output arrays &
            &have diffenrent shape'
       call terminate
    end if

    ! optical conductivity
    delt=0.d0
    if(oc(1).eq.oc(2)) delt=1.d0
    sigma(:)=aimag(eps(:))*w(:)/(4.d0*pi)
    sigma(:)=sigma(:)+zi*( -(dble(eps(:))-delt)*w(:)/(4.d0*pi) )

  end subroutine gensigma

end module m_gensigma
