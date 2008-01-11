subroutine calcupdatevectors(n,iunconverged,P,w,r,evalfv,evecfv,phi) 
  use modmain, only:nstfv,nmatmax,zzero,zone
  use diisinterfaces
  implicit none
  integer ,intent (in)::n , iunconverged
  complex(8),intent(in)::P(nmatmax,nmatmax)
  complex(8),intent(in)::r(n,nstfv),evecfv(n,nstfv)
  complex(8),intent(out)::phi(n,nstfv)
  real(8), intent(in)::w(nmatmax),evalfv(nstfv)
  complex(8)::z,v(n,nstfv)
  real(8)::nrm
  integer m,i,j
  m=nstfv
  call zgemm('C','N',n,m,n,zone,P,nmatmax,r,n,zzero,v,n)

  do j=1,m
     do i=1,n
        z=cmplx (w(i)-evalfv(j),0.0)
        if(abs(z).gt.1e-5)then  
           v(i,j)=-v(i,j)/z
        else
           v(i,j)=zzero
        endif
     end do

  end do

 
 do i=1,m
 call zcopy(n,evecfv(1,i),1,phi(1,i),1)
end do

  call zgemm('N','N',n,m,n,zone,P,nmatmax,v,n,zone,phi,n)
 do i=1,m
    ! call zaxpy(n,zone,phi(1,i),1,evecfv(1,i),1)
    call zcopy(n,phi(1,i),1,evecfv(1,i),1)

 end do

end subroutine calcupdatevectors
