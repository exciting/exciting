
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevecsv(vpl,evecsv)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer isym,lspn,ik,ist,i,j,k
integer recl,nstsv_
real(8) vkl_(3),det,th,t1,t2
real(8) s(3,3),sc(3,3),v(3)
complex(8) s2(2,2),zt1,zt2
! external functions
real(8) r3taxi
external r3taxi
! find the k-point number
call findkpt(vpl,isym,ik)
! index to global spin rotation in lattice point group
lspn=lspnsymc(isym)
! find the record length
inquire(iolength=recl) vkl_,nstsv_,evecsv
!$OMP CRITICAL
open(70,file=trim(scrpath)//'EVECSV'//trim(filext),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik) vkl_,nstsv_,evecsv
close(70)
!$OMP END CRITICAL
if (r3taxi(vkl(1,ik),vkl_).gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECSV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getevecsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EVECSV.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if symmetry element is the identity return
if (lspn.eq.1) return
! if eigenvectors are spin-unpolarised return
if (.not.spinpol) return
! spin rotation matrix
s(:,:)=dble(symlat(:,:,lspn))
! convert symmetry matrix to Cartesian coordinates
call r3mm(s,ainv,sc)
call r3mm(avec,sc,sc)
! determine the axis and angle of rotation for the symmetry matrix
call rotaxang(epslat,sc,det,v,th)
! determine the SU(2) representation of the symmetry matrix
t1=cos(th/2.d0)
t2=sin(th/2.d0)
s2(1,1)=t1
s2(1,2)=0.d0
s2(2,1)=0.d0
s2(2,2)=t1
do k=1,3
  zt1=-zi*t2*v(k)
  do i=1,2
    do j=1,2
      s2(i,j)=s2(i,j)+zt1*sigmat(i,j,k)
    end do
  end do
end do
! apply SU(2) symmetry matrix to second-variational states
do i=1,nstsv
  do ist=1,nstfv
    zt1=evecsv(ist,i)
    zt2=evecsv(ist+nstfv,i)
    evecsv(ist,i)=s2(1,1)*zt1+s2(1,2)*zt2
    evecsv(ist+nstfv,i)=s2(2,1)*zt1+s2(2,2)*zt2
  end do
end do
return
end subroutine

