
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: sumrule
! !INTERFACE:
subroutine sumrule(dynq)
! !INPUT/OUTPUT PARAMETERS:
!   dynq : dynamical matrices on q-point set (in,real(3*natmtot,3*natmtot,nqpt))
! !DESCRIPTION:
!   Applies the same correction to all the dynamical matrices such that the
!   matrix for ${\bf q}=0$ satisfies the acoustic sum rule. In other words, the
!   matrices are updated with
!   $$ D_{ij}^{\bf q}\rightarrow D_{ij}^{\bf q}-\sum_{k=1}^3 \omega_k^0
!    v_{k;i}^0 v_{k;j}^0 $$
!   for all ${\bf q}$, where $\omega_k^0$ is the $k$th eigenvalue of the
!   ${\bf q}=0$ dynamical matrix and $v_{k;i}^0$ the $i$th component of its
!   eigenvector. The eigenvalues are assumed to be arranged in ascending order.
!   This ensures that the ${\bf q}=0$ dynamical matrix has 3 zero eigenvalues,
!   which the uncorrected matrix may not have due to the finite
!   exchange-correlation grid.
!
! !REVISION HISTORY:
!   Created May 2005 (JKD)
!EOP
!BOC
use modmain
implicit none
! arguments
complex(8), intent(inout) :: dynq(3*natmtot,3*natmtot,nqpt)
! local variables
integer n,iq,iq0,i,j,k
integer lwork,info
real(8) t1
! allocatable arrays
real(8), allocatable :: w(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
complex(8), allocatable :: ev(:,:)
n=3*natmtot
allocate(w(n))
allocate(rwork(3*n))
lwork=2*n
allocate(work(lwork))
allocate(ev(n,n))
! find the Gamma-point phonon
do iq0=1,nqpt
  t1=sqrt(vql(1,iq0)**2+vql(2,iq0)**2+vql(3,iq0)**2)
  if (t1.lt.epslat) goto 10
end do
write(*,*)
write(*,'("Error(sumrule): no zero length q-vector")')
write(*,*)
stop
10 continue
! compute the eigenvalues and vectors of the q=0 dynamical matrix
do i=1,n
  do j=i,n
    ev(i,j)=0.5d0*(dynq(i,j,iq0)+conjg(dynq(j,i,iq0)))
  end do
end do
call zheev('V','U',n,ev,n,w,work,lwork,rwork,info)
! subtract outer products of 3 lowest eigenvectors for q=0 from all the
! dynamical matrices
do iq=1,nqpt
  do i=1,n
    do j=1,n
      do k=1,3
        dynq(i,j,iq)=dynq(i,j,iq)-w(k)*ev(i,k)*conjg(ev(j,k))
      end do
    end do
  end do
end do
deallocate(w,rwork,work,ev)
return
end subroutine
!EOC

