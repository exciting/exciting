

! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine getevecsv(vpl, evecsv)
use modmain
use modinput
use modmpi
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: evecsv(nstsv, nstsv)
! local variables
integer::isym, lspn, ik, ist, i
integer::recl, nstsv_
real(8)::vkl_(3), det, v(3), th, t1
complex(8) su2(2, 2), zt1, zt2
character(256) ::filetag
character(256), external:: outfilenamestring
!<chm>
logical::exist
integer::koffset
!</chm>
#ifdef XS
  ! added feature to access arrays for only a subset of bands
  complex(8), allocatable :: evecsv_(:, :)
#endif
! find the k-point number
call findkpt(vpl, isym, ik)
! index to global spin rotation in lattice point group
lspn=lspnsymc(isym)
! find the record length
#ifdef XS
inquire(iolength=recl) vkl_, nstsv_
#endif
#ifndef XS
inquire(iolength=recl) vkl_, nstsv_, evecsv
#endif
filetag=trim(filetag_evecsv)
!$OMP CRITICAL
do i=1, 100
inquire(file=outfilenamestring(filetag, ik), exist=exist)
 if (exist) then
open(70, file = outfilenamestring(filetag, ik), action = 'READ', &
 form = 'UNFORMATTED', access = 'DIRECT', recl = recl)
exit
else 
call system('sync')
write(*, *) "Waiting for other process to write"//":getevecsv:"// &
     trim(outfilenamestring(filetag, ik))
call sleep(5)
endif
enddo
 if (splittfile) then
 koffset=ik-firstk(procofk(ik))+1
 else
 koffset =ik
 endif
#ifdef XS
  read(70, rec=1) vkl_, nstsv_
  close(70)
  if (nstsv.gt.nstsv_) then
     write(*, *)
     write(*, '("Error(getevecsv): invalid nstsv for k-point ", I8)') ik
     write(*, '(" current    : ", I8)') nstsv
     write(*, '(" EVECSV.OUT : ", I8)') nstsv_
     write(*, '(" file	     : ", a	 )') trim(outfilenamestring(filetag, ik))
     write(*, *)
     stop
  end if
  allocate(evecsv_(nstsv_, nstsv_))
  inquire(iolength=recl) vkl_, nstsv_, evecsv_
  open(70, file = outfilenamestring(filetag, ik), action = 'READ', &
       form = 'UNFORMATTED', access = 'DIRECT', recl = recl)
  read(70, rec=koffset) vkl_, nstsv_, evecsv_
  ! retreive subset
  evecsv(:, :)=evecsv_(:nstsv, :nstsv)
  deallocate(evecsv_)
#endif
#ifndef XS
read(70, rec=koffset) vkl_, nstsv_, evecsv
#endif
close(70)
!$OMP END CRITICAL
t1=abs(vkl(1, ik)-vkl_(1))+abs(vkl(2, ik)-vkl_(2))+abs(vkl(3, ik)-vkl_(3))
if (t1.gt.input%structure%epslat) then
  write(*, *)
  write(*, '("Error(getevecsv): differing vectors for k-point ", I8)') ik
  write(*, '(" current	  : ", 3G18.10)') vkl(:, ik)
  write(*, '(" EVECSV.OUT : ", 3G18.10)') vkl_
  write(*, '(" file	  : ", a      )') trim(outfilenamestring(filetag, ik))
  write(*, *)
  stop
end if
#ifndef XS
if (nstsv.ne.nstsv_) then
  write(*, *)
  write(*, '("Error(getevecsv): differing nstsv for k-point ", I8)') ik
  write(*, '(" current	  : ", I8)') nstsv
  write(*, '(" EVECSV.OUT : ", I8)') nstsv_
  write(*, '(" file	  : ", a      )') trim(outfilenamestring(filetag, ik))
  write(*, *)
  stop
end if
#endif
! if symmetry element is the identity return
if (lspn.eq.1) return
! if eigenvectors are spin-unpolarised return
if (.not.associated(input%groundstate%spin)) return
! find the SU(2) representation of the spin rotation matrix
call rotaxang(input%structure%epslat, symlatc(:, :, lspn), det, v, th)
call axangsu2(v, th, su2)
! apply SU(2) symmetry matrix to second-variational states
do i=1, nstsv
  do ist=1, nstfv
    zt1=evecsv(ist, i)
    zt2=evecsv(ist+nstfv, i)
    evecsv(ist, i)=su2(1, 1)*zt1+su2(1, 2)*zt2
    evecsv(ist+nstfv, i)=su2(2, 1)*zt1+su2(2, 2)*zt2
  end do
end do
return
end subroutine

module m_getevecsvr 
  implicit none
contains


subroutine getevecsvr(fname, isti, istf, vpl, evecsv)
    use modmain
    implicit none
    ! arguments
    character(*), intent(in) :: fname
    integer, intent(in) :: isti, istf
    real(8), intent(in) :: vpl(3)
    complex(8), intent(out) :: evecsv(:, :)
    ! local variables
    integer :: err
    complex(8), allocatable :: evecsvt(:, :)
    character(256) :: str
    ! check correct shapes
    err=0
    if ((isti.lt.1).or.(istf.gt.nstsv).or.(istf.le.isti)) then
       write(*, *)
       write(*, '("Error(getevecsvr): inconsistent limits for bands:")')
       write(*, '(" band limits  : ", 2i6)') isti, istf
       write(*, '(" maximum value: ", i6)') nstsv
       write(*, *)
       err=err+1
    end if
    if ((size(evecsv, 1).ne.size(evecsv, 2)).or. &
	 (size(evecsv, 1).ne.(istf - isti + 1))) then
       write(*, *)
       write(*, '("Error(getevecsvr): output array does not match for bands:")')
       write(*, '(" band limits 	     : ", 2i6)') isti, istf
       write(*, '(" requested number of bands: ", i6)') istf-isti+1
       write(*, '(" array size		     : ", i6)') size(evecsv, 1)
       write(*, *)
       err=err+1
    end if
    if (err.ne.0) stop
    allocate(evecsvt(nstsv, nstsv))
    filetag_evecsv=trim(fname)
    str=trim(filext)
    filext=''
    call getevecsv(vpl, evecsvt)
    filetag_evecsv='EVECSV'
    filext=trim(str)
    evecsv(:, :)=evecsvt(isti:istf, isti:istf)
    deallocate(evecsvt)
  end subroutine getevecsvr
end module m_getevecsvr
