
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqnit(nmatp,ngp,igpig,vpl,vgpl,vgpc,apwalm,evalfv,evecfv)
use modmain
implicit none
! arguments
integer, intent(in) :: nmatp
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer is,ia,it,i
integer ist,jst
real(8) ts1,ts0
real(8) t1
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: h(:)
complex(8), allocatable :: o(:,:)
! external functions
complex(8) zdotc
external zdotc
call timesec(ts0)
allocate(o(nmatp,nstfv))
if ((iscl.ge.2).or.(task.eq.1).or.(task.eq.3)) then
! read in the eigenvalues/vectors from file
  call getevalfv(vpl,evalfv)
  call getevecfv(vpl,vgpl,evecfv)
else
! initialise the eigenvectors to canonical basis vectors
  evecfv(:,:)=0.d0
  do ist=1,nstfv
    evecfv(ist,ist)=1.d0
  end do
end if
! start iteration loop
do it=1,nseqit
! begin parallel loop over states
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(h,is,ia,t1,i)
!$OMP DO
  do ist=1,nstfv
    allocate(h(nmatp))
! operate with H and O on the current vector
    h(:)=0.d0
    o(:,ist)=0.d0
    do is=1,nspecies
      do ia=1,natoms(is)
        call hmlaa(.true.,is,ia,ngp,apwalm,evecfv(:,ist),h)
        call hmlalo(.true.,is,ia,ngp,apwalm,evecfv(:,ist),h)
        call hmllolo(.true.,is,ia,ngp,evecfv(:,ist),h)
        call olpaa(.true.,is,ia,ngp,apwalm,evecfv(:,ist),o(:,ist))
        call olpalo(.true.,is,ia,ngp,apwalm,evecfv(:,ist),o(:,ist))
        call olplolo(.true.,is,ia,ngp,evecfv(:,ist),o(:,ist))
      end do
    end do
    call hmlistl(.true.,ngp,igpig,vgpc,evecfv(:,ist),h)
    call olpistl(.true.,ngp,igpig,evecfv(:,ist),o(:,ist))
! normalise
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      do i=1,nmatp
        evecfv(i,ist)=t1*evecfv(i,ist)
        h(i)=t1*h(i)
        o(i,ist)=t1*o(i,ist)
      end do
    end if
! estimate the eigenvalue
    evalfv(ist)=dble(zdotc(nmatp,evecfv(:,ist),1,h,1))
! subtract the gradient of the Rayleigh quotient from the eigenvector
    t1=evalfv(ist)
    do i=1,nmatp
      evecfv(i,ist)=evecfv(i,ist)-tauseq*(h(i)-t1*o(i,ist))
    end do
! normalise
    o(:,ist)=0.d0
    do is=1,nspecies
      do ia=1,natoms(is)
        call olpaa(.true.,is,ia,ngp,apwalm,evecfv(:,ist),o(:,ist))
        call olpalo(.true.,is,ia,ngp,apwalm,evecfv(:,ist),o(:,ist))
        call olplolo(.true.,is,ia,ngp,evecfv(:,ist),o(:,ist))
      end do
    end do
    call olpistl(.true.,ngp,igpig,evecfv(:,ist),o(:,ist))
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      do i=1,nmatp
        evecfv(i,ist)=t1*evecfv(i,ist)
        o(i,ist)=t1*o(i,ist)
      end do
    end if
    deallocate(h)
! end parallel loop over states
  end do
!$OMP END DO
!$OMP END PARALLEL
! perform Gram-Schmidt orthonormalisation
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jst,zt1,t1,i)
!$OMP DO ORDERED
  do ist=1,nstfv
!$OMP ORDERED
    do jst=1,ist-1
      zt1=-zdotc(nmatp,evecfv(:,jst),1,o(:,ist),1)
      call zaxpy(nmatp,zt1,evecfv(:,jst),1,evecfv(:,ist),1)
      call zaxpy(nmatp,zt1,o(:,jst),1,o(:,ist),1)
    end do
!$OMP END ORDERED
! normalise
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      do i=1,nmatp
        evecfv(i,ist)=t1*evecfv(i,ist)
        o(i,ist)=t1*o(i,ist)
      end do
    end if
  end do
!$OMP END DO
!$OMP END PARALLEL
! end iteration loop
end do
deallocate(o)
call timesec(ts1)
!$OMP CRITICAL
timefv=timefv+ts1-ts0
!$OMP END CRITICAL
return
end subroutine

