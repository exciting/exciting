subroutine remove_converged(evecmap,iunconverged,rnorms,n,r,h,s,eigenvector,eigenvalue,trialvecs)
  use sclcontroll,only:diismax,epsresid,maxdiisspace
  use modmain,only:nstfv
  implicit none
  integer, intent(in)::n
  integer, intent(inout) ::evecmap(nstfv),iunconverged
  complex(8),intent(inout)::r(n,nstfv),h(n,nstfv,diismax),s(n,nstfv,diismax)
  complex(8),intent(inout)::eigenvector(n,nstfv),trialvecs(n,nstfv,diismax)
  real(8),intent(inout)::eigenvalue(nstfv,diismax),rnorms(nstfv)
  integer i,skipp,idiis,oldindex,newindex
  skipp=0
  do i=1,nstfv
     if (evecmap(i).gt.0) then
        if(rnorms(evecmap(i)).lt.epsresid)then
           evecmap(i)=0	
           skipp=skipp+1
        else
           if (skipp.gt.0) then
              oldindex=evecmap(i)
              newindex=evecmap(i)-skipp	
              call zcopy(n,eigenvector(1,oldindex),1,eigenvector(1,newindex),1)
              call zcopy(n,r(1,oldindex),1,r(1,newindex),1)
             
              rnorms(newindex)=rnorms(oldindex)
              do idiis=1,maxdiisspace
               eigenvalue(newindex,idiis)=eigenvalue(oldindex,idiis)
                 call zcopy(n,h(1,oldindex,idiis),1,h(1,newindex,idiis),1)
                 call zcopy(n,s(1,oldindex,idiis),1,s(1,newindex,idiis),1)
                 call zcopy(n,trialvecs(1,oldindex,idiis),1,trialvecs(1,newindex,idiis),1)
              end do
              evecmap(i)=newindex
           endif
        endif
     endif
  end do
  iunconverged=0
  do i=1,nstfv
     if (evecmap(i).gt.0) then
        iunconverged=iunconverged+1
     endif
  end do
 ! write(*,*)iunconverged,"map",evecmap
end subroutine remove_converged

