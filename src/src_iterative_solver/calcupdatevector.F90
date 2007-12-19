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
  integer m,i
  m=nstfv
  !call zgemm('C','N',n,m,n,zone,P,nmatmax,r,n,zzero,v,n)
  
  do i=1,m
   call zgemv('C',n,n,zone,P,nmatmax,r(1,i),1,zzero,v(1,i),1)
     if(abs(w(i)-evalfv(i)).lt.1e-2)then
        z= cmplx(1.0/(w(i)-evalfv(i)),0)
        call zscal(n,z,v(1,i),1)
        write(*,*)"hier z",z
       call zgemv('N',n,n,zone,P,nmatmax,v(1,i),1,zzero,v(1,i),1) 
     else
      v(:,i)=0
     endif
  end do
  do i=1,m
  	call zcopy(n,evecfv(1,i),1,phi(1,i),1)
  end do
   write(771,*),phi(:,2)

  call zgemm('N','N',n,m,n,zone,P,nmatmax,v,n,zone,phi,n)
 
   write(772,*),phi(:,2) 
end subroutine calcupdatevectors
