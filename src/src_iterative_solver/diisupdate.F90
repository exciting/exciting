subroutine   diisupdate(idiis,iunconverged,n,h,s,trialvec,evalfv ,evecfv)
  use modmain,only: nstfv,nmatmax,zone,zzero
  use diisinterfaces
  implicit none
  integer ,intent(in)::idiis,iunconverged,n
  complex(8),intent(in)::h(n,nstfv,idiis),s(n,nstfv,idiis),trialvec(n,nstfv,idiis)
  real(8), intent(in):: evalfv(nstfv)
  complex(8),intent(out)::evecfv(nmatmax,nstfv)


  complex(8) p(n,idiis)
  real(8)::nrm
  integer::i,j,ir,is
  complex(8):: Pmatrix(idiis,idiis), Qmatrix(idiis,idiis),c(idiis)
   complex(8)::z
    evecfv(:,:)=0.0
  do i=1,iunconverged 

     do j=1,idiis
        call zcopy(n,h(1,i,j),1,p(1,j),1)
        z=cmplx(-evalfv(i),0)
        call zaxpy(n,z,s(1,i,j),1,p(1,j),1)
     end do
     call zgemm('C','N',idiis,idiis,n,zone,p,n,p,n,&
         zzero,Pmatrix,idiis)
     do ir=1,idiis
	do is=1,idiis
           Qmatrix(is,ir)=zdotc(n,trialvec(1,i,is),1,s(1,i,ir),1)
	enddo
     enddo
     call solvediis(idiis,Pmatrix,Qmatrix,c)
    
     do ir=1,idiis 
        call zaxpy(n,-c(ir),trialvec(1,i,ir),1,evecfv(1,i),1)
        nrm=sqrt( dble( zdotc( n,evecfv(1,i),1,evecfv(1,i),1 ) ) )
        z=cmplx(1.0/nrm,0)
        call zscal(n,z,evecfv(1,i),1)
     end do

  end do
end subroutine diisupdate
