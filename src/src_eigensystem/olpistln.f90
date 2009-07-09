


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olpistl
! !INTERFACE:


subroutine olpistln(overlap, ngp, igpig)
! !USES:
use modmain
use modfvsystem
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
type(HermiteanMatrix), intent(inout)::overlap
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)


! local variables
integer::i, j, k, iv(3), ig
complex(8) zt1

! calculate the matrix elements
!$omp parallel default(shared) & 
!$omp  private(iv,ig,i,j)
!$omp do
  do j=1, ngp
    !k=((j-1)*j)/2
    do i=1, j
      !k=k+1
      iv(:)=ivg(:, igpig(i))-ivg(:, igpig(j))
      ig=ivgig(iv(1), iv(2), iv(3))
      if ((ig.gt.0).and.(ig.le.ngvec)) then
	  call Hermiteanmatrix_indexedupdate(overlap, j, i, cfunig(ig))
      end if
    end do
  end do
!$omp end do 
!$omp end parallel

return
end subroutine
!EOC
