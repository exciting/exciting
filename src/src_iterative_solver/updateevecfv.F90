subroutine  updateevecfv(n,m,da,nmatmax,evecfv,evalfv,evecp,evalp)
	implicit none
integer, intent(in)::n,m,nmatmax
complex(8), intent(in)::da(n,m),evecp(2*m,m)
complex(8),intent(inout)::evecfv(nmatmax,m)
real(8),intent (in)::evalp(2*m)
real(8),intent (out)::evalfv(m)

!local vars
complex(8)::basis(n,2*m)
do i=1,2*m 
	if(i.le.m) then
	call zcopy(n,evecfv(:,i),1,basis(:,i),1)
	else
	call zcopy(n,da(:,i-m),1,basis(:,i),1)
	endif
end do
evecfv=0
do i=1,m
	do j=1,2*m
		call zaxpy(n,evecp(j,i),basis(:,j),1,evecfv(:,i),1)
		!evecfv(:,i)=evecfv(:,i)+basis(:,j)*evecp(j,i)
	end do
end do
do i=1,m
	evalfv(i)=evalp(i)
end do

end subroutine
