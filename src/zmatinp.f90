
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: zmatinp
! !INTERFACE:
subroutine zmatinp(tapp,n,alpha,x,y,v,a)
! !INPUT/OUTPUT PARAMETERS:
!   tapp  : .true. if the matrix is to be applied to the input vector v,
!           .false. if the full matrix is to be calculated (in,logical)
!   n     : length of vectors (in,integer)
!   alpha : complex constant (in,complex)
!   x     : first input vector (in,complex(n))
!   y     : second input vector (in,complex(n))
!   v     : input vector to which matrix is applied if tapp is .true., otherwise
!           not referenced (in,complex(n))
!   a     : matrix applied to v if tapp is .true., otherwise the full matrix in
!           packed form (inout,complex(n+(n-1)*n/2))
! !DESCRIPTION:
!   Performs the rank-2 operation
!   $$ A_{ij}\rightarrow\alpha{\bf x}_i^*{\bf y}_j+\alpha^*{\bf y}_i^*{\bf x}_j
!    +A_{ij}, $$
!   where $A$ is stored in packed form. This is similar to the {\tt BLAS}
!   routine {\tt zhpr2}, except that here a matrix of inner products is formed
!   instead of an outer product of vectors. If {\tt tapp} is {\tt .true.} then
!   the matrix is applied to an input vector, rather than calculated explicitly.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: n
complex(8), intent(in) :: alpha
complex(8), intent(in) :: x(n)
complex(8), intent(in) :: y(n)
complex(8), intent(in) :: v(n)
complex(8), intent(inout) :: a(*)
! local variables
integer i,j,k
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-12
complex(8) zt1,zt2
if (tapp) then
! apply the matrix
  zt1=y(1)*v(1)
  zt2=x(1)*v(1)
  do j=2,n
    zt1=zt1+y(j)*v(j)
    zt2=zt2+x(j)*v(j)
  end do
  zt1=conjg(alpha*zt1)
  zt2=alpha*conjg(zt2)
  if ((abs(dble(zt1)).gt.eps).or.(abs(aimag(zt1)).gt.eps).or. &
   (abs(dble(zt2)).gt.eps).or.(abs(aimag(zt2)).gt.eps)) then
    do i=1,n
      a(i)=a(i)+conjg(zt1*x(i)+zt2*y(i))
    end do
  end if
else
! calculate the matrix elements
  k=0
  do j=1,n
    if ((abs(dble(x(j))).gt.eps).or.(abs(aimag(x(j))).gt.eps).or. &
     (abs(dble(y(j))).gt.eps).or.(abs(aimag(y(j))).gt.eps)) then
      zt1=conjg(alpha*y(j))
      zt2=alpha*conjg(x(j))
      do i=1,j-1
        k=k+1
        a(k)=a(k)+conjg(zt1*x(i)+zt2*y(i))
      end do
      k=k+1
      a(k)=dble(a(k))+2.d0*dble(zt1*x(j))
    else
      k=k+j
    end if
  end do
end if
return
end subroutine
!EOC
