
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevalfv(vpl,evalfv)
use modmain
use modmpi
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: evalfv(nstfv,nspnfv)
! local variables
integer isym,ik,koffset,i
logical exist
integer recl,nstfv_,nspnfv_
real(8) vkl_(3),t1

character(256) ::filetag
character(256), external:: outfilenamestring
#ifdef XS
  ! added feature to access arrays for only a subset of bands
  real(8), allocatable :: evalfv_(:,:)
#endif
! find the k-point number
call findkpt(vpl,isym,ik)
! find the record length

#ifdef XS
inquire(iolength=recl) vkl_,nstfv_,nspnfv_
#endif
#ifndef XS
inquire(iolength=recl) vkl_,nstfv_,nspnfv_,evalfv
#endif

filetag='EVALFV'
do i=1,100
inquire(file=outfilenamestring(filetag,ik),exist=exist)
 if (exist) then
open(70,file=outfilenamestring(filetag,ik),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
exit
else 
call system('sync')
write(*,*) "Waiting for other process to write"//":getevalfv:"// &
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
  read(70,rec=1) vkl_,nstfv_,nspnfv_
  close(70)
  if (nstfv.gt.nstfv_) then
     write(*,*)
     write(*,'("Error(getevalfv): invalid nstfv for k-point ",I8)') ik
     write(*,'(" current    : ",I8)') nstfv
     write(*,'(" EVALFV.OUT : ",I8)') nstfv_
     write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
     write(*,*)
     stop
  end if
  allocate(evalfv_(nstfv_,nspnfv_))
  inquire(iolength=recl) vkl_,nstfv_,nspnfv_,evalfv_
  open(70,file=outfilenamestring(filetag,ik),action='READ', &
       form='UNFORMATTED',access='DIRECT',recl=recl)
  read(70,rec=koffset) vkl_,nstfv_,nspnfv_,evalfv_
  ! retreive subset
  evalfv(:,:)=evalfv_(:nstfv,:)
  deallocate(evalfv_)
#endif

#ifndef XS
read(70,rec=koffset) vkl_,nstfv_,nspnfv_,evalfv
#endif
close(70)



t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevalfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVALFV.OUT : ",3G18.10)') vkl_
  write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
  write(*,*)
  stop
end if
#ifndef XS
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstfv
  write(*,'(" EVALFV.OUT : ",I8)') nstfv_
  write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
  write(*,*)
  stop
end if
#endif
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspnfv
  write(*,'(" EVALFV.OUT : ",I8)') nspnfv_
  write(*,'(" file       : ",a      )') trim(outfilenamestring(filetag,ik))
  write(*,*)
  stop
end if
return
end subroutine

module m_getevalfvr
  implicit none
contains
  subroutine getevalfvr(isti,istf,vpl,evalfv)
    use modmain
    implicit none
    ! arguments
    integer, intent(in) :: isti,istf
    real(8), intent(in) :: vpl(3)
    real(8), intent(out) :: evalfv(:,:)
    ! local variables
    integer :: err
    real(8), allocatable :: evalfvt(:,:)
    ! check correct shapes
    err=0
    if ((isti.lt.1).or.(istf.gt.nstfv).or.(istf.le.isti)) then
       write(*,*)
       write(*,'("Error(getevalfvr): inconsistent limits for bands:")')
       write(*,'(" band limits  : ",2i6)') isti,istf
       write(*,'(" maximum value: ",i6)') nstfv
       write(*,*)
       err=err+1
    end if
    if (size(evalfv,1).ne.(istf-isti+1)) then
       write(*,*)
       write(*,'("Error(getevalfvr): output array does not match for bands:")')
       write(*,'(" band limits              : ",2i6)') isti,istf
       write(*,'(" requested number of bands: ",i6)') istf-isti+1
       write(*,'(" array size               : ",i6)') size(evalfv,1)
       write(*,*)
       err=err+1
    end if
    if (size(evalfv,2).ne.nspnfv) then
       write(*,*)
       write(*,'("Error(getevalfvr): output array does not match for &
            &nspnfv:")')
       write(*,'(" nspnfv    : ",i6)') nspnfv
       write(*,'(" array size: ",i6)') size(evalfv,2)
       write(*,*)
       err=err+1       
    end if
    if (err.ne.0) stop
    allocate(evalfvt(nstfv,nspnfv))
    call getevalfv(vpl,evalfvt)
    evalfv(:,:)=evalfvt(isti:istf,:)
    deallocate(evalfvt)
  end subroutine getevalfvr
end module m_getevalfvr
