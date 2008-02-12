
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genwiq2xs
! !INTERFACE:
subroutine genwiq2xs(flag,iq,igq1,igq2,clwt)
  ! !USES:
  use modmain
  use modxs
  use m_genfilname
  use m_getunit
  ! !DESCRIPTION:
  !   Effective integrals of Coulomb interaction.
  !
  ! !REVISION HISTORY:
  !   Created August 2008 (SAG)
  !EOP
  !BOC
  implicit none
  ! arguments
  integer, intent(in) :: flag,iq,igq1,igq2
  real(8), intent(out) :: clwt
  ! local variables
  integer, parameter :: np=4
  integer, parameter :: ns0=10,nss=20
  integer ns,i1,i2,i3,i,ip
  real(8) d(3),dv,sum2,t1,t2
  real(8) v0(3),v1(3),v2(3),v3(3)
  real(8) xa(np),ya(np),c(np)
  integer :: un,n,j1,j2
  real(8) :: v0l(3),v1l(3),v01(3),v02(3),v11(3)
  real(8) :: v12(3),v21(3),v22(3),v31(3),v32(3)
  real(8) :: sum3,t3,qsz,cpu0,cpu1
  character(256) :: fname
  ! external functions
  real(8) polynom
  external polynom
  call cpu_time(cpu0)

!!$  ! map the q-vector into the first Brillouin zone
!!$  t1=1.d8
!!$  v0(:)=0.d0
!!$  do i1=-1,1
!!$     do i2=-1,1
!!$        do i3=-1,1
!!$           v1(:)=vqc(:,iq)+dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2) &
!!$                +dble(i3)*bvec(:,3)
!!$           v1l(:)=vql(:,iq)+(/dble(i1),dble(i2),dble(i3)/)
!!$           t2=v1(1)**2+v1(2)**2+v1(3)**2
!!$           if (t2.lt.t1) then
!!$              t1=t2
!!$              v0(:)=v1(:)
!!$              v0l(:)=v1l(:)
!!$           end if
!!$        end do
!!$     end do
!!$  end do

  ! do not map to the first Brillouin zone
  v0(:)=vqc(:,iq)
  v0l(:)=vql(:,iq)

  ! integrate out singularites and/or improve accuracy in summations
  if (flag.eq.0) then
     ! modified Coulomb potential in reciprocal space
     clwt=fourpi/(gqc(igq1,iq)*gqc(igq2,iq))
  else if (flag.eq.1) then
     ! radius of sphere with same volume than subcell
     qsz=(6*pi**2/(omega*nqpt))**(1.d0/3.d0)
     if ((igq1.ne.1).or.(igq2.ne.1)) then
        ! integrate out 1/q singularity by spherical Volume
        clwt=(qsz**2*omega/pi )/gqc(igq2,iq)
     else
        ! integrate out 1/q^2 singularity by spherical Volume
        clwt=2*qsz*omega/pi
     end if
  else if (flag.eq.2) then
     j1=igqig(igq1,iq)
     j2=igqig(igq2,iq)
     v01(:)=v0(:)+dble(ivg(:,j1))
     v02(:)=v0(:)+dble(ivg(:,j2))
     ! loop over different subdivisions
     ns=ns0
     do ip=1,np
        ! subdivision vectors in lattice coordinates
        do i=1,3
           d(i)=1.d0/(dble(ngridq(i)*2*ns))
        end do
        ! smallest volume element
        dv=((twopi**3)/omega)*d(1)*d(2)*d(3)
        ! compute the integral of 1/q^2
        sum2=0.d0
        do i1=-ns,ns-1
           t1=dble(i1)*d(1)
           v11(:)=v01(:)+t1*bvec(:,1)
           v12(:)=v02(:)+t1*bvec(:,1)
           do i2=-ns,ns-1
              t1=dble(i2)*d(2)
              v21(:)=v11(:)+t1*bvec(:,2)
              v22(:)=v12(:)+t1*bvec(:,2)
              do i3=-ns,ns-1
                 t1=dble(i3)*d(3)
                 v31(:)=v21(:)+t1*bvec(:,3)
                 v32(:)=v22(:)+t1*bvec(:,3)
                 t3=sqrt(sum(v31**2)*sum(v32**2))
                 if (t3.gt.1.d-14) then
                    sum3=sum3+1.d0/t3
                 end if
              end do
           end do
        end do
        sum2=sum2*dv
        sum3=sum3*dv
        xa(ip)=dv**(1.d0/3.d0)
        ya(ip)=sum3 !+-
        ! increment number of subdivisions
        ns=ns+nss
     end do
     ! extrapolate the volume element to zero with a polynomial
     clwt=polynom(0,np,xa,ya,c,0.d0)*fourpi
  end if
  call cpu_time(cpu1)
  t1=cpu1-cpu0
!!$  if (flag.ne.0) write(*,'("Time in genwiq2xs (seconds):",f12.2)') t1
end subroutine genwiq2xs
!EOC

