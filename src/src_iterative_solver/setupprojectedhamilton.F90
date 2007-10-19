subroutine setupprojectedhamilton(n,m,h,nmatmax,evecfv,da,hprojected,oprojected)
	implicit none
integer,intent(in)::n,m,nmatmax
complex(8), intent(in)::evecfv(nmatmax,m),h(n*(n+1)/2),da(n,m)
complex(8), intent(out)::hprojected(2*m*(2*m+1)/2),oprojected(2*m*(2*m+1)/2)
complex(8):: basis(n,2*m),vec(n)
complex(8) ::zdotc
external zdotc
integer::i,j,pi
do i=1,2*m 
	if(i.le.m) then
		call zcopy(n,evecfv(:,i),1,basis(:,i),1)
	else
		call zcopy(n,da(:,i-m),1,basis(:,i),1)
	endif
end do
pi=1
do i=1,2*m
	do j=1,i
		vec=0.0
		call zhpmv('U',n,(1.d0,1.d0),h,basis(1,j), 1, (0.d0,0.d0), vec(1), 1)
		hprojected(pi)=zdotc(n ,basis(1,i),1,vec(1),1)
		oprojected(pi)=zdotc(n ,basis(1,i),1,basis(1,j),1)
		if(i.eq.j) then
			hprojected(pi)=dble(hprojected(pi))
		endif
		pi=pi+1
	end do
end do

end subroutine
