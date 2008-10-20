
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
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
!   Effective integrals of Coulomb interaction. Based on routine {\tt genwiq2}.
!
! !REVISION HISTORY:
!   Created February 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: flag,iq,igq1,igq2
  real(8), intent(out) :: clwt
  ! local variables
  integer, parameter :: ns0=10,nss=20
  integer :: np,ns,i1,i2,i3,i,ip,nrbox
  real(8) :: d(3),dv,sum2,t1,blim(2),blen,vllim(3),ran(3),ranl(3),omegabox
  real(8) :: t2,qsz,cpu0,cpu1
  real(8), allocatable :: xa(:),ya(:),c(:)
  real(8) :: v01(3),v02(3),v11(3),v12(3),v21(3),v22(3),v31(3),v32(3)
  ! external functions
  real(8) polynom
  external polynom
  call cpu_time(cpu0)
  ! determine G+q-vectors
  v01(:)=vgqc(:,igq1,iq)
  v02(:)=vgqc(:,igq2,iq)
  ! integrate out singularites and/or improve accuracy in summations
  select case(flag)
  case(0)
     ! modified Coulomb potential in reciprocal space
     clwt=fourpi/(gqc(igq1,iq)*gqc(igq2,iq))
  case(1)
     ! radius of sphere with same volume than subcell
     qsz=(6*pi**2/(omega*nqpt))**(1.d0/3.d0)
     if ((igq1.ne.1).or.(igq2.ne.1)) then
        ! integrate out 1/q singularity by spherical Volume
        clwt=(qsz**2*omega*nqpt/pi)/gqc(igq2,iq)
     else if ((igq1.eq.1).and.(igq2.eq.1)) then
        ! integrate out 1/q^2 singularity by spherical Volume
        clwt=2*qsz*omega*nqpt/pi
     else
        write(*,*)
        write(*,'("Error(genwiq2xs): analytic method chosen for regular case")')
        write(*,*)
        call terminate
     end if
  case(2)
     np=2
     ! higher order extrapolation for 1/q^2 term
     if ((igq1.eq.1).and.(igq2.eq.1)) np=3
     allocate(xa(np),ya(np),c(np))
     ! loop over different subdivisions
     ns=ns0
     do ip=1,np
        ! subdivision vectors in lattice coordinates
        do i=1,3
           d(i)=1.d0/(dble(ngridq(i)*2*ns))
        end do
        ! smallest volume element (we drop the (2pi)^3/omega factor!)
        dv=d(1)*d(2)*d(3)
        ! compute the integral of 1/(|q+G||q+Gp|) over the small
        ! fraction of the Brillouin zone centered at q 
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
                 t2=sqrt(sum(v31**2)*sum(v32**2))
                 if (t2.gt.1.d-14) then
                    sum2=sum2+1.d0/t2
                 end if
              end do
           end do
        end do
        sum2=sum2*dv
        xa(ip)=dv**(1.d0/3.d0)
        ya(ip)=sum2
        ! increment number of subdivisions
        ns=ns+nss
     end do
     ! extrapolate the volume element to zero with a polynomial
     clwt=polynom(0,np,xa,ya,c,0.d0)*fourpi*nqpt
     deallocate(xa,ya,c)
  case(3)
     ! RIM method inspired by Andrea Marini
     ! find maximum extension of small Brillouin zone
     blim(:)=0.d0
     do i1=-1,1,2
        v11(:)=dble(i1)*bvec(:,1)/dble(ngridq(1)*2)
        do i2=-1,1,2
           v21(:)=v11(:)+dble(i2)*bvec(:,2)/dble(ngridq(2)*2)
    	   do i3=-1,1,2
	      v31(:)=v21(:)+dble(i3)*bvec(:,3)/dble(ngridq(3)*2)
	      do i=1,3
	         ! lower limit for box
                 if (v31(i).lt.blim(1)) blim(1)=v31(i)
		 ! upper limit for box
		 if (v31(i).gt.blim(2)) blim(2)=v31(i)
              end do
	   end do
        end do
     end do
     ! box length
     blen=blim(2)-blim(1)
     ! limits of sBZ in lattice      
     vllim(:)=1.d0/dble(2*ngridq(:))
     ! needs high value (above 1e6) to converge - inferior to method Nr. 2
     nrbox=1000000
     t1=0.d0
     omegabox=blen**3
     do i=1,10000
        call random_number(t2)
     end do
     do i=1,nrbox
        call random_number(ran)
	! map random number to box
	ran(:)=ran(:)*blen
	ran(:)=ran(:)+blim(1)
	! check if random vector is in sBZ
	ranl=matmul(binv,ran)
	if (all(ranl.gt.-vllim).and.(all(ranl.lt.vllim)).and.&
             (sum(abs(ran)).gt.1.d-14)) then
	   t2=sum((v01+ran)**2)*sum((v02+ran)**2)
	   t1=t1+1.d0/sqrt(t2)	   
	end if
     end do
     clwt=t1*(omegabox/nrbox)*fourpi*nqpt*omega/(twopi)**3
  end select
  call cpu_time(cpu1)
  t1=cpu1-cpu0
!!$  if (flag.ne.0) write(*,'("Time in genwiq2xs (seconds):",3i6,f12.2)') iq, &
!!$       igq1,igq2,t1
end subroutine genwiq2xs
!EOC

