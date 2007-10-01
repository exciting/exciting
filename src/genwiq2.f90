
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genwiq2
! !INTERFACE:
subroutine genwiq2
! !USES:
use modmain
! !DESCRIPTION:
!   The Fock matrix elements
!   $$ V^{\rm NL}_{ij{\bf k}}\equiv\sum_{l{\bf k'}}\int
!    \frac{\phi^*_{i{\bf k}}({\bf r})\phi_{l{\bf k}'}({\bf r})
!    \phi^*_{l{\bf k}'}({\bf r}')\phi_{j{\bf k}}({\bf r}')}{|{\bf r}-{\bf r'}|}
!    \,d{\bf r}\,d{\bf r'} $$
!   contain a divergent term in the sum over ${\bf k}'$ which behaves as
!   $1/q^2$, where ${\bf q}\equiv{\bf k}-{\bf k}'$ is in the first Brillouin
!   zone. The resulting convergence with respect to the number of discrete
!   $q$-points is very slow. This routine computes the weights
!   \begin{align}\label{genwiq2_1}
!    w_{{\bf q}_i}\equiv\int_{V_i}\frac{1}{q^2}\,d{\bf q}\;,
!   \end{align}
!   where the integral is over the small parallelepiped centered on ${\bf q}_i$,
!   so that integrals over the first Brillouin zone of the form
!   $$ I=\int_{\rm BZ}\frac{f({\bf q})}{q^2}\,d{\bf q}\;, $$
!   can be approximated by the sum
!   $$ I\approx\sum_i w_{{\bf q}_i}f({\bf q}_i) $$
!   which converges rapidly with respect to the number of $q$-points for smooth
!   functions $f$. The integral in (\ref{genwiq2_1}) is determined by evaluating
!   it numerically on increasingly finer grids and extrapolating to the
!   continuum. Agreement with Mathematica to at least 10 significant figures.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD,SS)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: np=5
integer, parameter :: ns0=10,nss=20
integer ns,iq,i1,i2,i3,i,ip
real(8) d(3),dv,sum,t1,t2
real(8) v0(3),v1(3),v2(3),v3(3)
real(8) xa(np),ya(np),c(np)
! external functions
real(8) polynom
external polynom
! allocate global wiq2 array
if (allocated(wiq2)) deallocate(wiq2)
allocate(wiq2(ngridq(1)*ngridq(2)*ngridq(3)))
! if system is a molecule wiq2 should be zero
if (molecule) then
  wiq2(:)=0
  return
end if
! begin loop over q-points
do iq=1,nqpt
! map the q-vector into the first Brillouin zone
  t1=1.d8
  v0(:)=0.d0
  do i1=-1,1
    do i2=-1,1
      do i3=-1,1
        v1(:)=vqc(:,iq)+dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2) &
         +dble(i3)*bvec(:,3)
        t2=v1(1)**2+v1(2)**2+v1(3)**2
        if (t2.lt.t1) then
          t1=t2
          v0(:)=v1(:)
        end if
      end do
    end do
  end do
! loop over different subdivisions
  ns=ns0
  do ip=1,np
! subdivision vectors in lattice coordinates
    do i=1,3
      d(i)=1.d0/(dble(ngridk(i)*2*ns))
    end do
! smallest volume element
    dv=((twopi**3)/omega)*d(1)*d(2)*d(3)
! compute the integral of 1/q^2
    sum=0.d0
    do i1=-ns,ns-1
      t1=dble(i1)*d(1)
      v1(:)=v0(:)+t1*bvec(:,1)
      do i2=-ns,ns-1
        t1=dble(i2)*d(2)
        v2(:)=v1(:)+t1*bvec(:,2)
        do i3=-ns,ns-1
          t1=dble(i3)*d(3)
          v3(:)=v2(:)+t1*bvec(:,3)
          t2=v3(1)**2+v3(2)**2+v3(3)**2
          if (t2.gt.1.d-14) then
            sum=sum+1.d0/t2
          end if
        end do
      end do
    end do
    sum=sum*dv
    xa(ip)=dv**(1.d0/3.d0)
    ya(ip)=sum
! increment number of subdivisions
    ns=ns+nss
  end do
! extrapolate the volume element to zero with a polynomial
  wiq2(iq)=polynom(0,np,xa,ya,c,0.d0)
end do
return
end subroutine
!EOC

