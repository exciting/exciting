subroutine   diisupdate(idiis,iunconverged,n,h,s,trialvec,evalfv ,evecfv,infodiisupdate)
  use modmain,only: nstfv,zone,zzero
  use diisinterfaces
  use sclcontroll,only:recalculate_preconditioner
  use sclcontroll,only:diismax
  implicit none
  integer ,intent(in)::idiis,iunconverged,n
  complex(8),intent(in)::h(n,nstfv, diismax),s(n,nstfv, diismax),trialvec(n,nstfv, diismax)
  real(8), intent(in):: evalfv(nstfv)
  complex(8),intent(out)::evecfv(n,nstfv)
  integer, intent(out)::infodiisupdate
  


  complex(8) p(n,idiis)
  real(8)::nrm
  integer::i,j,ir,is
real(8):: Pmatrix(idiis+1,idiis+1), Qmatrix(idiis+1,idiis+1),c(idiis+1),residnorm2
   complex(8)::z
   
  do i=1,iunconverged 
  !calculate residuals
     do j=1,idiis
        call zcopy(n,h(1,i,j),1,p(1,j),1)
        z=cmplx(-evalfv(i),0)
        call zaxpy(n,z,s(1,i,j),1,p(1,j),1)
     end do
    residnorm2=dble(zdotc(n,p(1,idiis),1,p(1,idiis),1))
     do ir=1,idiis
		do is=1,idiis
			Pmatrix(is,ir)=dble(zdotc(n,p(1,is),1,p(1,ir),1))/residnorm2
			!if (dble(Pmatrix(is,ir)).lt.1.e-4) then
			!write(889,*)"ir,is,p(1,i,ir)",ir,is,p(1,is)
			!endif
	 	enddo
     enddo
     do ir=1,idiis
		do is=1,idiis
			Qmatrix(is,ir)=dble(zdotc(n,trialvec(1,i,is),1,s(1,i,ir),1))
			if (dble(Qmatrix(is,ir)).lt.1.e-4) then
			write(888,*)"ir,is,trialvec(1,i,is),s(1,i,ir)",ir,is,trialvec(:,i,is),"s\n\n",s(:,i,ir)
			endif
			
		enddo
     enddo
     call solvediislin(idiis,Pmatrix,Qmatrix,c)
     if  (recalculate_preconditioner .eqv. .true.) then
     infodiisupdate=1
     exit
     endif
    write(*,*) "c",c
    evecfv(:,i)=0.0
     do ir=1,idiis
     z=cmplx(c(ir),0.0)
        call zaxpy(n, z,trialvec(1,i,ir),1,evecfv(1,i),1)
     end do
  end do
   infodiisupdate=0
end subroutine diisupdate
