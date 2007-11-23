
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
integer isym,lspn,ik,ist,i
integer recl,nstsv_
real(8) vkl_(3),det
complex(8) su2(2,2),zt1,zt2
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
! find the SU(2) representation of the spin rotation matrix
call rotsu2(symlatc(1,1,lspn),det,su2)
! apply SU(2) symmetry matrix to second-variational states
do i=1,nstsv
  do ist=1,nstfv
    zt1=evecsv(ist,i)
    zt2=evecsv(ist+nstfv,i)
    evecsv(ist,i)=su2(1,1)*zt1+su2(1,2)*zt2
    evecsv(ist+nstfv,i)=su2(2,1)*zt1+su2(2,2)*zt2
  end do
end do
return
end subroutine

