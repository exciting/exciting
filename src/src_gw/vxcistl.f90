
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: vxcistl
! !INTERFACE:
subroutine vxcistl(tapp,ngp,igpig,vgpc,v,h)
! !USES:
use modmain
use modgw
! !INPUT/OUTPUT PARAMETERS:
!   tapp  : .true. if the Hamiltonian is to be applied to the input vector,
!           .false. if the full matrix is to be calculated (in,logical)
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   v     : input vector to which H is applied if tapp is .true., otherwise
!           not referenced (in,complex(nmatmax))
!   h     : H applied to v if tapp is .true., otherwise it is the Hamiltonian
!           matrix in packed form (inout,complex(npmatmax))
! !DESCRIPTION:
!   Computes the interstitial contribution to the Hamiltonian matrix for the APW
!   basis functions. The Hamiltonian is given by
!   $$ H^{\rm I}({\bf G+k,G'+k})=\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}({\bf G-G'})+V^{\sigma}({\bf G-G'}), $$
!   where $V^{\sigma}$ is the effective interstitial potential for spin $\sigma$
!   and $\tilde{\Theta}$ is the characteristic function. See routine
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: v(nmatmax)
complex(8), intent(inout) :: h(*)
! local variables
integer i,j,k,ig,iv(3)
real(8) t1
complex(8) zt1
if (tapp) then
! apply the Hamiltonian operator to v
  do i=1,ngp
    do j=i,ngp
      iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
      iv(:)=modulo(iv(:)-intgv(:,1),ngrid(:))+intgv(:,1)
      ig=ivgig(iv(1),iv(2),iv(3))
      if ((ig.gt.0).and.(ig.le.ngvec)) then
        zt1=vxcig(ig)
        h(i)=h(i)+zt1*v(j)
        if (i.ne.j) h(j)=h(j)+conjg(zt1)*v(i)
      end if
    end do
  end do
else
! calculate the matrix elements
  do j=1,ngp
    k=((j-1)*j)/2
    do i=1,j
      k=k+1
      iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
      iv(:)=modulo(iv(:)-intgv(:,1),ngrid(:))+intgv(:,1)
      ig=ivgig(iv(1),iv(2),iv(3))
      if ((ig.gt.0).and.(ig.le.ngvec)) then
        h(k)=h(k)+vxcig(ig)
      end if
    end do
  end do
end if
return
end subroutine
!EOC
