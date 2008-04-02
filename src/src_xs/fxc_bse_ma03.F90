
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_fxc_bse_ma03
  implicit none

contains

  subroutine fxc_bse_ma03(msiz,sw,iw,fxc)
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
    integer, intent(in) :: iw
    complex(8), intent(out) :: fxc(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_bse_ma03'
    character(256) :: filnam
    complex(8) :: zt1
    integer :: sh(2),ig,igp,un,a1,a2,a3,recl

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
    
    do ig=1,msiz
       do igp=1,msiz
          read(un,rec=iw) fxc
       end do
    end do

    if (.not.sw) then
       zt1=fxc(1,1)
       fxc(:,:)=zzero
       fxc(1,1)=zt1
    end if


    ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! apply scaling for debugging purposes
    fxc = conjg(fxc)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!write(*,*) 'reading fxc-bse:',iw,fxc
!fxc(1,1)=dble(fxc(1,1))
!fxc(1,1)=-0.2/fourpi
!fxc(1,1)=(-8.191520659266058E-002,-3.652876264792654E-003)
!fxc(1,1)=(-8.191520659266058E-002,0.d0)



    close(un)

  end subroutine fxc_bse_ma03

end module m_fxc_bse_ma03
