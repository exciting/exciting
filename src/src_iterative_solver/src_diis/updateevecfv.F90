



subroutine  updateevecfv(n, m, da, nmatmax, evecfv, evalfv, evecp, evalp)
	implicit none
integer, intent(in)::n, m, nmatmax
complex(8), intent(in)::da(n, m), evecp(2*m, m)
complex(8), intent(inout)::evecfv(nmatmax, m)
real(8), intent (in)::evalp(2*m)
real(8), intent (out)::evalfv(m)
integer::i, j

!local vars
complex(8)::basis(n, 2*m)
do i=1, 2*m 
	if(i.le.m) then
	call zcopy(n, evecfv(1, i), 1, basis(1, i), 1)
	else
	call zcopy(n, da(1, i-m), 1, basis(1, i), 1)
	endif
end do
evecfv(:, :)=0.d0
do i=1, m
	do j=1, 2*m
		call zaxpy(n, evecp(j, i), basis(1, j), 1, evecfv(1, i), 1)
		!evecfv(:,i)=evecfv(:,i)+basis(:,j)*evecp(j,i)
	end do
end do

end subroutine
