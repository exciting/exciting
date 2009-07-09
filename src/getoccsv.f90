


! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine getoccsv(vpl, occsvp)
  use modmain
use modinput
  use modmpi
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3)
  real(8), intent(out) :: occsvp(nstsv)
  ! local variables
  logical::exist
  integer::isym, ik, koffset, i
  integer::recl, nstsv_
  real(8)::vkl_(3), t1
  character(256) ::filetag
  character(256), external:: outfilenamestring
#ifdef XS
  ! added feature to access arrays for only a subset of bands
  real(8), allocatable :: occsv_(:)
#endif
  ! find the k-point number
  call findkpt(vpl, isym, ik)
  ! find the record length

#ifdef XS
  inquire(iolength=recl) vkl_, nstsv_
#endif
#ifndef XS
  inquire(iolength=recl) vkl_, nstsv_, occsvp
#endif
  filetag=trim(filetag_occsv)
 do i=1, 100
inquire(file=outfilenamestring(filetag, ik), exist=exist)
 if (exist) then
open(70, file = outfilenamestring(filetag, ik), action = 'READ', &
 form = 'UNFORMATTED', access = 'DIRECT', recl = recl)
exit
else 
call system('sync')
write( * , *) "Waiting for other process to write"//":getoccsv:"// &
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
     write(*, '("Error(getoccsv): invalid nstsv for k-point ", I8)') ik
     write(*, '(" current    : ", I8)') nstsv
     write(*, '(" OCCSV.OUT  : ", I8)') nstsv_
     write(*, '(" file	     : ", a	 )') trim(outfilenamestring(filetag, ik))
     write(*, *)
     stop
  end if
  allocate(occsv_(nstsv_))
  inquire(iolength=recl) vkl_, nstsv_, occsv_
  open(70, file = outfilenamestring(filetag, ik), action = 'READ', &
       form = 'UNFORMATTED', access = 'DIRECT', recl = recl)
  read(70, rec=koffset) vkl_, nstsv_, occsv_
  ! retreive subset
  occsvp(:)=occsv_(:nstsv)
  deallocate(occsv_)
#endif
#ifndef XS
read(70, rec=koffset) vkl_, nstsv_, occsvp
#endif
  close(70)

t1=abs(vkl(1, ik)-vkl_(1))+abs(vkl(2, ik)-vkl_(2))+abs(vkl(3, ik)-vkl_(3))
  if (t1.gt.input%structure%epslat) then
     write(*, *)
     write(*, '("Error(getoccsv): differing vectors for k-point ", I8)') ik
     write(*, '(" current    : ", 3G18.10)') vkl(:, ik)
     write(*, '(" OCCSV.OUT  : ", 3G18.10)') vkl_
     write(*, '(" file	     : ", a	 )') trim(outfilenamestring(filetag, ik))
     write(*, *)
     stop
  end if
#ifndef XS
  if (nstsv.ne.nstsv_) then
     write(*, *)
     write(*, '("Error(getoccsv): differing nstsv for k-point ", I8)') ik
     write(*, '(" current    : ", I8)') nstsv
     write(*, '(" OCCSV.OUT  : ", I8)') nstsv_
     write(*, '(" file	     : ", a	 )') trim(outfilenamestring(filetag, ik))
     write(*, *)
     stop
  end if
#endif
  return
end subroutine getoccsv

module m_getoccsvr
  implicit none
contains


subroutine getoccsvr(fname, isti, istf, vpl, occsvp)
    use modmain
    implicit none
    ! arguments
    character(*), intent(in) :: fname
    integer, intent(in) :: isti, istf
    real(8), intent(in) :: vpl(3)
    real(8), intent(out) :: occsvp(:)
    ! local variables
    integer :: err
    real(8), allocatable :: occsvt(:)
    character(256) :: str
    ! check correct shapes
    err=0
    if ((isti.lt.1).or.(istf.gt.nstsv).or.(istf.le.isti)) then
       write(*, *)
       write(*, '("Error(getoccsvr): inconsistent limits for bands:")')
       write(*, '(" band limits  : ", 2i6)') isti, istf
       write(*, '(" maximum value: ", i6)') nstsv
       write(*, *)
       err=err+1
    end if
    if (size(occsvp, 1).ne.(istf-isti+1)) then
       write(*, *)
       write(*, '("Error(getoccsvr): output array does not match for bands:")')
       write(*, '(" band limits 	     : ", 2i6)') isti, istf
       write(*, '(" requested number of bands: ", i6)') istf-isti+1
       write(*, '(" array size		     : ", i6)') size(occsvp, 1)
       write(*, *)
       err=err+1
    end if
    if (err.ne.0) stop
    allocate(occsvt(nstsv))
    filetag_occsv=trim(fname)
    str=trim(filext)
    filext=''
    call getoccsv(vpl, occsvt)
    filetag_occsv='OCCSV'
    filext=trim(str)
    occsvp(:)=occsvt(isti:istf)
    deallocate(occsvt)
  end subroutine getoccsvr
end module m_getoccsvr
