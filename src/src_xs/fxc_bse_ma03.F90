
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_fxc_bse_ma03
  implicit none

contains

  subroutine fxc_bse_ma03(msiz,sw,w,fxc)
    !
    ! BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
    !
    !
    use modmain
    use modmpi
    use modxs
    use invert
    use m_tdgauntgen
    use m_findgntn0
    use m_writegqpts
    use m_genfilname
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: msiz
    ! true if all G-components of fxc are to be considered
    logical, intent(in) :: sw
    complex(8), intent(in) :: w
    complex(8), intent(out) :: fxc(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_bse_ma03'
    complex(8) :: zt1
    integer :: sh(2),ig

    sh=shape(fxc)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout,'(a,2i9,a,i9,a)') 'Error('//trim(thisnam)//'): size of &
            &fxc is to small (required)', sh, '(', msiz, ')'
       call terminate
    end if


    ! *** read kernel from file ***
    ! KERNXC_BSE.OUT
    fxc(:,:)=zzero

  end subroutine fxc_bse_ma03

end module m_fxc_bse_ma03
