subroutine   diisupdate(idiis,iunconverged,n,h,s,trialvec,evalfv ,evecfv)
use modmain,only: nstfv,nmatmax
implicit none
integer ,intent(in)::idiis,iunconverged,n
complex(8),intent(in)::h(n,nstfv,idiis),s(n,nstfv,idiis),trialvec(n,nstfv,idiis)
real(8), intent(in):: evalfv(nstfv)
complex(8),intent(out)::evecfv(nmatmax,nstfv)
complex(8) zdotc
real(8) dlamch
external zdotc,dlamch
complex p(n,idiis)
real(8)::rnorms(idiis)
integer::i,j,ir,is
complex(8):: Pmatrix(idiis,idiis), Qmatrix(idiis,idiis),c(idiis)
do i=1,iunconverged 
	do j=1,idiis
	call residualvectors(n,1,h(:,i,j),s(:,i,j),&
	trialvec(:,i,j),p(:,j),rnorms)
	end do
	!Pmatrix :gemm(pT.p)
	call zgemm('C','N',n,idiis,n,idiis,complex(1,0),p,n,p,n,complex(0,0),Pmatrix,idiis)
	do ir=1,idiis
	do is=1,idiis
	Qmatrix(ir,is)=zdotc(n,trialvec(1,i,ir),1,s(1,i,is),1)
	enddo
	enddo	
	call solvediis(idiis,Pmatrix,Qmatrix,c)
	do
	!evecfv(:,i): gemv(subspacevectors(:,:,i) .alpha)
	end do

end do
end subroutine
