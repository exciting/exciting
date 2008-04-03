subroutine exactupdatevectors(n,iunconverged,hamilton,overlap,r,&
rhizvalue,eigenvector,trialvecs)  
!calculate update equation with linsolver   

!solvefor dA:  dA=(H-e*S)\R

! dA 	Update step to zero residual
! H 	Hamilton 
! S 	Overlap
! e 	Rhitz Value
! R 	Residual

!trialvecs=eigenvector+dA
use modfvsystem
use modmain,only:zone
integer, intent(in):: n,iunconverged
type(HermiteanMatrix),intent(in)::hamilton,overlap
complex(8),intent(in)::r(n,iunconverged),eigenvector(n,iunconverged)
real(8),intent(in)::rhizvalue(iunconverged)
complex(8),intent(out)::trialvecs(n,iunconverged)
complex(8)::dA(n,iunconverged),sigma
type(HermiteanMatrix):: HES
call zcopy(n*iunconverged,eigenvector, 1 ,trialvecs,1)
call zcopy(n*iunconverged,r,1, dA,1)
call newmatrix(HES,ispacked(hamilton),n)
Do i=1,iunconverged
	call HermiteanMatrixcopy(hamilton,HES)
	sigma=dcmplx(-rhizvalue(i),0)
	call HermiteanMatrixAXPY(sigma,overlap,HES)
	call HermiteanmatrixLU(HES)
	call Hermiteanmatrixlinsolve(hes,dA(:,i))
	call zaxpy(n,zone,dA(1,i),1,trialvecs(1,i),1)
end do
call deletematrix(HES)
end subroutine