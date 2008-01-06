
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ematqkgir2
  implicit none
contains

  subroutine ematqkgir2(iq,ik,igq)
    use modmain
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,igq
    ! local variables
    character(*), parameter :: thisnam = 'ematqkgir2'
    integer :: ikq,ig,ig1,ig2,ig3,igk0,igk,iv(3),iv1(3),iv3(3)
    ! wrapping of BZ
    integer :: ivwrap(3)
    integer, allocatable :: aigk0(:),aigk(:)

    ikq=ikmapikq(ik,iq)
    allocate(aigk0(ngkmax0),aigk(ngkmax))

    ! positive wrapping G-vector
    ivwrap(:)=nint(vkl0(:,ik)+vql(:,iq)-vkl(:,ikq))

    ! precalculate for speedup
    aigk0(:)=igkig0(:,ik,1)
    aigk(:)=igkig(:,ikq,1)
    ig3=igqig(igq,iq)
    iv3(:)=ivg(:,ig3)
    do igk0 = 1, ngk0(ik,1)
       ig1=aigk0(igk0)
       iv1(:)=ivg(:,ig1)+iv3(:)
       do igk = 1, ngk(ikq,1)
          ig2=aigk(igk)
          ! wrapping of k+q vector included
          iv(:)=iv1(:)-(ivg(:,ig2) - ivwrap(:))
          ig = ivgig(iv(1),iv(2),iv(3))
          xihir(igk0,igk) = cfunig(ig)
       end do
    end do

    deallocate(aigk0,aigk)

  end subroutine ematqkgir2
  
end module m_ematqkgir2
