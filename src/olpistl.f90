
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olpistl
! !INTERFACE:
subroutine olpistl(tapp,ngp,igpig,v,o)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   v     : input vector to which O is applied if tapp is .true., otherwise
!           not referenced (in,complex(nmatmax))
!   o     : O applied to v if tapp is .true., otherwise it is the overlap
!         : matrix in packed form (inout,complex(npmatmax))
! !DESCRIPTION:
!   Computes the interstitial contribution to the overlap matrix for the APW
!   basis functions. The overlap is given by
!   $$ O^{\rm I}({\bf G+k,G'+k})=\tilde{\Theta}({\bf G-G'}), $$
!   where $\tilde{\Theta}$ is the characteristic function. See routine
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
complex(8), intent(in) :: v(nmatmax)
complex(8), intent(inout) :: o(*)
! local variables
integer i,j,k,iv(3),ig
complex(8) zt1
if (tapp) then
! apply the overlap operator to v
  do i=1,ngp
    do j=i,ngp
      iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
      ig=ivgig(iv(1),iv(2),iv(3))
      if ((ig.gt.0).and.(ig.le.ngvec)) then
        zt1=cfunig(ig)
        o(i)=o(i)+zt1*v(j)
        if (i.ne.j) o(j)=o(j)+conjg(zt1)*v(i)
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
      ig=ivgig(iv(1),iv(2),iv(3))
      if ((ig.gt.0).and.(ig.le.ngvec)) then
        o(k)=o(k)+cfunig(ig)
      end if
    end do
  end do
end if
return
end subroutine
!EOC

