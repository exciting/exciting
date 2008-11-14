subroutine minenergy(sigma)
use modmain
implicit none
complex(8),intent(out)::sigma

real(8)::lorbearray(maxlorbord*maxlorb*maxspecies)
integer::en,i,j,k,m
integer,external::idamax
real(8),parameter::one=1
 m=1
  lorbearray=0
 en =maxlorbord*maxlorb*maxspecies
 do i=1,maxlorbord
 do j=1,maxlorb
 do k=1,maxspecies
 if (abs (lorbe0(i,j,k)).lt.20)  lorbearray(m)=lorbe0(i,j,k)
 m=m+1
 enddo
 enddo
 enddo
call dscal(en,-one,lorbearray,1)
sigma=dcmplx( min( -abs(lorbearray(idamax(en,lorbearray,1)))  ,-1.d0))
write(*,*)sigma
end subroutine
