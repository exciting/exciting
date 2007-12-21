subroutine calcupdatevectors(n,iunconverged,P,w,r,evalfv,phi) 
  use modmain, only:nstfv,nmatmax,zzero,zone
  use diisinterfaces
  implicit none
  integer ,intent (in)::n , iunconverged
  complex(8),intent(in)::P(nmatmax,nmatmax),r(n,nstfv)
  complex(8),intent(out)::phi(n,nstfv)
  real(8), intent(in)::w(nmatmax),evalfv(nstfv)
  complex(8):: v(n,nstfv)
  complex(8)::z
  integer m,i
  m=nstfv
  call zgemm('C','N',n,m,n,zone,P,nmatmax,r,n,zzero,v,n)

  do i=1,m
     if(abs(w(i)-evalfv(i)).lt.1e-1)then
   		z= cmplx(1.0/(w(i)-evalfv(i)),0)
        call zscal(n,z,v,1)
        write(*,*)"hier"
     else
        v(:,i)=0
     endif
  end do

  call zgemm('N','N',n,m,n,zone,P,nmatmax,v,n,zzero,phi,n)

end subroutine calcupdatevectors
