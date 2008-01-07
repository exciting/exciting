
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
    real(8) :: vkql(3)
    integer :: isym,ikq,ig,ig1,ig2,ig3,igk0,igk,iv(3),iv1(3),iv3(3),ivu(3)
    integer, allocatable :: aigk0(:),aigk(:)
    ! k-point index for equivalent point to k+q
    vkql(:)=vkl(:,ik)+vql(:,iq)
    call findkpt(vkql,isym,ikq)
    allocate(aigk0(ngkmax),aigk(ngkmax))
    ! positive wrapping G-vector
    ivu(:)=nint(vkl(:,ik)+vql(:,iq)-vkl(:,ikq))
    ! precalculate for speedup
    aigk0(:)=igkig(:,ik,1)
    aigk(:)=igkig(:,ikq,1)
    ig3=igqig(igq,iq)
    iv3(:)=ivg(:,ig3)
    do igk0=1,ngk(ik,1)
       ig1=aigk0(igk0)
       iv1(:)=ivg(:,ig1)+iv3(:)
       do igk=1,ngk(ikq,1)
          ig2=aigk(igk)
          ! wrapping of k+q vector included
          iv(:)=iv1(:)-(ivg(:,ig2)-ivu(:))
          ig = ivgig(iv(1),iv(2),iv(3))
          xihir(igk0,igk)=cfunig(ig)
       end do
    end do
    deallocate(aigk0,aigk)
  end subroutine ematqkgir2
  
end module m_ematqkgir2
