! function to compose filename for parallel execution
! REMARK: never call with stringconstat allways call by reference to filetag

character(256) function outfilenamestring(filetag,ik)
use modmpi, only:procs,lastk,firstk,procofk,splittfile
use modmain, only:scrpath,filext
implicit none
!external lastk,firstk

!character(256):: outfilenamestring
character(256), intent(in) :: filetag
integer,intent(in)::ik
character(256):: tmp,tmp2,krange,scrpathtmp
 krange=''
 tmp=''
 tmp2=''
 scrpathtmp=''
 outfilenamestring=''
#ifdef MPI

 if ((procs.gt.1).and.splittfile) then
     write(tmp,'(I5)')firstk(procofk(ik))
     write(tmp2,'(I5)')lastk(procofk(ik))
     krange=trim(adjustl(tmp))//'-'//trim(adjustl(tmp2))
     scrpathtmp=scrpath
  endif
#endif
outfilenamestring=trim(scrpathtmp)//trim(filetag)//trim(krange)//trim(filext)   
end function outfilenamestring
