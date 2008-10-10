
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
integer is,ia,it,i,m
integer ist,jst,info
integer iwork(5*2),ifail(2)
real(8) ts0,ts1
real(8) vl,vu,w(2),t1
real(8) rwork(7*2)
complex(8) ap(3),bp(3),z(2,2)
complex(8) work(2,2),zt1,zt2
! allocatable arrays
complex(8), allocatable :: h(:)
complex(8), allocatable :: o(:,:)
complex(8), allocatable :: g(:)
complex(8), allocatable :: hg(:)
complex(8), allocatable :: og(:)
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
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(h,g,hg,og) &
!$OMP PRIVATE(is,ia,t1,ap,bp,vl,vu) &
!$OMP PRIVATE(m,w,z,work,rwork) &
!$OMP PRIVATE(iwork,ifail,info)
!$OMP DO
  do ist=1,nstfv
    allocate(h(nmatp),g(nmatp),hg(nmatp),og(nmatp))
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
! estimate the eigenvalue
    evalfv(ist)=dble(zdotc(nmatp,evecfv(:,ist),1,h,1))
! limit magnitude of eigenvalue
    if (abs(evalfv(ist)).gt.5.d0) evalfv(ist)=0.d0
! compute the gradient of the Rayleigh quotient
    t1=evalfv(ist)
    g(:)=h(:)-t1*o(:,ist)
! operate with H and O on the current gradient
    hg(:)=0.d0
    og(:)=0.d0
    do is=1,nspecies
      do ia=1,natoms(is)
        call hmlaa(.true.,is,ia,ngp,apwalm,g,hg)
        call hmlalo(.true.,is,ia,ngp,apwalm,g,hg)
        call hmllolo(.true.,is,ia,ngp,g,hg)
        call olpaa(.true.,is,ia,ngp,apwalm,g,og)
        call olpalo(.true.,is,ia,ngp,apwalm,g,og)
        call olplolo(.true.,is,ia,ngp,g,og)
      end do
    end do
    call hmlistl(.true.,ngp,igpig,vgpc,g,hg)
    call olpistl(.true.,ngp,igpig,g,og)
! solve the 2x2 generalised eigenvalue problem in the basis {g,evecfv}
    ap(1)=zdotc(nmatp,g,1,hg,1)
    ap(2)=zdotc(nmatp,g,1,h,1)
    ap(3)=evalfv(ist)
    bp(1)=zdotc(nmatp,g,1,og,1)
    bp(2)=zdotc(nmatp,g,1,o(:,ist),1)
    bp(3)=zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1)
    call zhpgvx(1,'V','I','U',2,ap,bp,vl,vu,1,1,-1.d0,m,w,z,2,work,rwork, &
     iwork,ifail,info)
    zt1=z(1,1)
    zt2=z(2,1)
    do i=1,nmatp
      evecfv(i,ist)=zt1*g(i)+zt2*evecfv(i,ist)
    end do
    deallocate(h,g,hg,og)
  end do
!$OMP END DO
!$OMP END PARALLEL
! perform Gram-Schmidt orthonormalisation
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,ia,t1,i,jst,zt1)
!$OMP DO
  do ist=1,nstfv
    do jst=1,ist-1
      zt1=-zdotc(nmatp,evecfv(:,jst),1,o(:,ist),1)
      call zaxpy(nmatp,zt1,evecfv(:,jst),1,evecfv(:,ist),1)
    end do
! normalise
    t1=abs(dble(zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1)))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      do i=1,nmatp
        evecfv(i,ist)=t1*evecfv(i,ist)
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

