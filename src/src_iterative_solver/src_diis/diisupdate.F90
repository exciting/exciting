subroutine   diisupdate(idiis,iunconverged,n,h,s&
     ,trialvec,evalfv ,evecfv,infodiisupdate)
  use modmain,only: nstfv,zone,zzero
  use diisinterfaces
  use sclcontroll,only:recalculate_preconditioner
  use sclcontroll,only:diismax
  implicit none
  integer ,intent(in)::idiis,iunconverged,n
  complex(8),intent(in)::h(n,nstfv, diismax)
  complex(8),intent(in):: s(n,nstfv, diismax),trialvec(n,nstfv, diismax)
  real(8), intent(in):: evalfv(nstfv,diismax)
  complex(8),intent(out)::evecfv(n,nstfv)
  integer, intent(out)::infodiisupdate
  logical lin	


  complex(8) p(n,idiis)
  real(8)::nrm
  integer::i,j,ir,is,idx(idiis)
  real(8):: Pmatrix(idiis+1,idiis+1), Qmatrix(idiis+1,idiis+1)
  real(8)::   c(idiis+1),residnorm2
  complex(8)::z
  lin=.true.
  Pmatrix=0.0
  Qmatrix=0.0
  p=0
  c=0

  
  do i=1,iunconverged 
     !calculate residuals
     do j=1,idiis
        call zcopy(n,h(1,i,j),1,p(1,j),1)
        z=cmplx(-evalfv(i,j),0)
        call zaxpy(n,z,s(1,i,j),1,p(1,j),1)
     end do
     residnorm2=dble(zdotc(n,p(1,idiis),1,p(1,idiis),1))

     do ir=1,idiis
        do is=1,idiis
           Pmatrix(is,ir)=dble(zdotc(n,p(1,is),1,p(1,ir),1))/ residnorm2
           if (dble(Pmatrix(is,ir)).lt.1.e-4) then
           write(889,*)"ir,is,p(1,i,ir)",ir,is,p(1,is)
           endif
        enddo
     enddo
     if(.not.lin) then
        do ir=1,idiis
           do is=1,idiis
              Qmatrix(is,ir)=dble(zdotc(n,trialvec(1,i,is),1,s(1,i,ir),1))
              if (dble(Qmatrix(is,ir)).lt.1.e-4) then
                 ! write(*,*)"warning Qmatrix(is,ir)).lt.1.e-4 in diisupdate"
                 !   write(888,*)"ir,is,trialvec(1,i,is),s(1,i,ir)",ir,is,trialvec(:,i,is),"s\n\n",s(:,i,ir)
              endif

           enddo
        enddo
     endif
    ! if(i==1 )write(*,*)"Pmatrix",Pmatrix
     if(lin) then
        call solvediislin(idiis,Pmatrix,Qmatrix,c)
         call sortidx(idiis,-abs(c),idx)
        if(i==1 )write(*,*)"c:",c(idx),"sum",sum(c(idx))
        c=c/sum(c(idx))
     else
		call solvediis(idiis,Pmatrix,Qmatrix,c)
     endif
     if  (recalculate_preconditioner .eqv. .true.) then
        infodiisupdate=1
        exit
     endif
     ! write(*,*) "c",c
   
   
call zcopy(n,zzero,0,evecfv(:,i),1)
     do ir=1,idiis
        z=dcmplx(c(idx(ir)),0.0)
        call zaxpy(n, z,trialvec( 1,i, idx(ir) ) ,1,evecfv(1,i),1)
        if(i==1 )write (*,*)"trialnrm",i,dznrm2( n,  trialvec( 1,i, idx(ir)),1)
     end do
  end do
  call zcopy(n*iunconverged,evecfv,1,trialvec( 1,1, idiis ),1)
  infodiisupdate=0
end subroutine diisupdate
