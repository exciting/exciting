subroutine   diisupdate(idiis,iunconverged,n,h,s,trialvec,evalfv ,evecfv)
use modmain,only: nstfv,nmatmax
implicit none
integer ,intent(in)::idiis,iunconverged,n
complex(8),intent(in)::h(n,nstfv,idiis),s(n,nstfv,idiis),trialvec(n,nstfv,idiis)
real(8), intent(in):: evalfv(nstfv)
complex(8),intent(out)::evecfv(nmatmax,nstfv)
complex p(n,idiis)
real(8)::rnorms(idiis)
integer::i,j
do i=1,iunconverged 
do j=1,idiis
call residualvectors(n,1,h(:,j,i),s(:,j,i),&
trialvec(:,j,i),p(:,j),rnorms)
end do
!Pmatrix :gemmv(pT.p)
!Qmatrix :gemmv(subspacecvectors(:,:,i)T.s)
!alpha:solve (P,Q system)
!evecfv(:,i): gemmv(subspacevectors(:,:,i) .alpha)
end do
end subroutine
