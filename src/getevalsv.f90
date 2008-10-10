
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
  real(8) vkl_(3),t1
    character(256) ::filetag
  character(256), external:: outfilenamestring

#ifdef XS
  ! added feature to access arrays for only a subset of bands
  real(8), allocatable :: evalsv_(:)
#endif
  ! find the k-point number
  call findkpt(vpl,isym,ik)

  ! find the record length
#ifndef XS
  inquire(iolength=recl) vkl_,nstsv_
#endif
#ifdef XS
  inquire(iolength=recl) vkl_,nstsv_,evalsvp
#endif
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
#ifdef XS
  read(70,rec=1) vkl_,nstsv_
  close(70)
  if (nstsv.gt.nstsv_) then
     write(*,*)
     write(*,'("Error(getevalsv): invalid nstsv for k-point ",I8)') ik
     write(*,'(" current    : ",I8)') nstsv
     write(*,'(" EVALSV.OUT : ",I8)') nstsv_
     write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
     write(*,*)
     stop
  end if
  allocate(evalsv_(nstsv_))
  inquire(iolength=recl) vkl_,nstsv_,evalsv_
  open(70,file=outfilenamestring(filetag,ik),action='READ', &
       form='UNFORMATTED',access='DIRECT',recl=recl)
  read(70,rec=koffset) vkl_,nstsv_,evalsv_
  ! retreive subset
  evalsvp(:)=evalsv_(:nstsv)
  deallocate(evalsv_)
#endif
#ifndef XS
  read(70,rec=koffset) vkl_,nstsv_,evalsvp
#endif
close(70)

t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
  if (t1.gt.epslat) then
     write(*,*)
     write(*,'("Error(getevalsv): differing vectors for k-point ",I8)') ik
     write(*,'(" current    : ",3G18.10)') vkl(:,ik)
     write(*,'(" EVALSV.OUT : ",3G18.10)') vkl_
     write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
     write(*,*)
     stop
  end if
#ifndef XS
  if (nstsv.ne.nstsv_) then
     write(*,*)
     write(*,'("Error(getevalsv): differing nstsv for k-point ",I8)') ik
     write(*,'(" current    : ",I8)') nstsv
     write(*,'(" EVALSV.OUT : ",I8)') nstsv_
     write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
     write(*,*)
     stop
  end if
#endif
  return
end subroutine getevalsv

module m_getevalsvr
  implicit none
contains
  subroutine getevalsvr(isti,istf,vpl,evalsvp)
    use modmain
    implicit none
    ! arguments
    integer, intent(in) :: isti,istf
    real(8), intent(in) :: vpl(3)
    real(8), intent(out) :: evalsvp(:)
    ! local variables
    integer :: err
    real(8), allocatable :: evalsvt(:)
    ! check correct shapes
    err=0
    if ((isti.lt.1).or.(istf.gt.nstsv).or.(istf.le.isti)) then
       write(*,*)
       write(*,'("Error(getevalsvr): inconsistent limits for bands:")')
       write(*,'(" band limits  : ",2i6)') isti,istf
       write(*,'(" maximum value: ",i6)') nstsv
       write(*,*)
       err=err+1
    end if
    if (size(evalsvp,1).ne.(istf-isti+1)) then
       write(*,*)
       write(*,'("Error(getevalsvr): output array does not match for bands:")')
       write(*,'(" band limits              : ",2i6)') isti,istf
       write(*,'(" requested number of bands: ",i6)') istf-isti+1
       write(*,'(" array size               : ",i6)') size(evalsvp,1)
       write(*,*)
       err=err+1
    end if
    if (err.ne.0) stop
    allocate(evalsvt(nstsv))
    call getevalsv(vpl,evalsvt)
    evalsvp(:)=evalsvt(isti:istf)
    deallocate(evalsvt)
  end subroutine getevalsvr
end module m_getevalsvr
