
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevalsv(vpl,evalsvp)
  use modmain
  use modmpi
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3)
  real(8), intent(out) :: evalsvp(nstsv)
  ! local variables
  logical exist
  integer isym,ik,koffset,i
  integer recl,nstsv_
  real(8) vkl_(3)
  ! external functions
  real(8) r3taxi
    character(256) ::filetag
  character(256), external:: outfilenamestring

  external r3taxi
  ! find the k-point number
  call findkpt(vpl,isym,ik)

  ! find the record length
  inquire(iolength=recl) vkl_,nstsv_,evalsvp
  filetag='EVALSV'
 do i=1,100
inquire(file=outfilenamestring(filetag,ik),exist=exist)
 if (exist) then
open(70,file=outfilenamestring(filetag,ik),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
exit
else 
call system('sync')
write(*,*) "Waiting for other process to write"//":getevalsv:"// &
     trim(outfilenamestring(filetag,ik))
call sleep(5)
endif
enddo
  if (splittfile) then
     koffset=ik-firstk(procofk(ik))+1
  else
     koffset =ik
  endif
  read(70,rec=koffset) vkl_,nstsv_,evalsvp
close(70)

  if (r3taxi(vkl(1,ik),vkl_).gt.epslat) then
     write(*,*)
     write(*,'("Error(getevalsv): differing vectors for k-point ",I8)') ik
     write(*,'(" current    : ",3G18.10)') vkl(:,ik)
     write(*,'(" EVALSV.OUT : ",3G18.10)') vkl_
     write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
     write(*,*)
     stop
  end if
  if (nstsv.ne.nstsv_) then
     write(*,*)
     write(*,'("Error(getevalsv): differing nstsv for k-point ",I8)') ik
     write(*,'(" current    : ",I8)') nstsv
     write(*,'(" EVALSV.OUT : ",I8)') nstsv_
     write(*,*)
     stop
  end if
  return
end subroutine getevalsv

