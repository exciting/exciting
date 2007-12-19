subroutine calcupdatevectors(n,iunconverged,P,w,r,evalfv,evecfv,phi) 
  use modmain, only:nstfv,nmatmax,zzero,zone
  use diisinterfaces
  implicit none
  integer ,intent (in)::n , iunconverged
  complex(8),intent(in)::P(nmatmax,nmatmax)
  complex(8),intent(in)::r(n,nstfv),evecfv(nmatmax,nstfv)
  complex(8),intent(out)::phi(n,nstfv)
  real(8), intent(in)::w(nmatmax),evalfv(nstfv)
  complex(8)::z,v(n,nstfv)
  integer m,i,j
  m=nstfv
  call zgemm('C','N',n,m,n,zone,P,nmatmax,r,n,zzero,v,n)

  do j=1,m
     do i=1,n
        z= w(i)-evalfv(j)
        if(dble(z).lt.1e-2)then  
           v(i,j)=v(i,j)/z
        else
           v(i,j)=0
        endif
     end do

  end do

 
 do i=1,m
 call zcopy(n,evecfv(1,i),1,phi(1,i),1)
end do

  call zgemm('N','N',n,m,n,zone,P,nmatmax,v,n,zzero,phi,n)
 do i=1,m
     call zaxpy(n,zone,phi(1,i),1,evecfv(1,i),1)
 end do

end subroutine calcupdatevectors
