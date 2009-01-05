
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_fxc_bse_ma03
  implicit none

contains

!BOP
! !ROUTINE: fxc_bse_ma03
! !INTERFACE:
  subroutine fxc_bse_ma03(msiz,sw,iw,fxc)
! !USES:
    use modmain
    use modmpi
    use modxs
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
    integer, intent(in) :: msiz
    ! true if all G-components of fxc are to be considered
    logical, intent(in) :: sw
    integer, intent(in) :: iw
    complex(8), intent(out) :: fxc(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_bse_ma03'
    character(256) :: filnam
    complex(8) :: zt1
    integer :: sh(2),un,recl
    sh=shape(fxc)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout,'(a,2i9,a,i9,a)') 'Error('//trim(thisnam)//'): size of &
            &fxc is to small (required)', sh, '(', msiz, ')'
       call terminate
    end if
    ! filename for BSE-xc-kernel
    call getunit(un)
    ! filename for xc-kernel
    call genfilname(basename='FXC_BSE',asc=.false.,bzsampl=bzsampl,&
         acont=acont,nar=.not.aresdf,iqmt=1,filnam=filnam)
    inquire(iolength=recl) fxc(:,:)
    open(un,file=trim(filnam),form='unformatted',action='read', &
         status='old',access='direct',recl=recl)
    read(un,rec=iw) fxc
    if (.not.sw) then
       zt1=fxc(1,1)
       fxc(:,:)=zzero
       fxc(1,1)=zt1
    end if
    close(un)
  end subroutine fxc_bse_ma03
!EOC

end module m_fxc_bse_ma03
