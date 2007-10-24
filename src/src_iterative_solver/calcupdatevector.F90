subroutine calcupdatevector(n,residual,HminuseS,da)
	implicit none
	integer, intent(in)::n
    complex(8),intent(in)::residual(n)
	complex(8),intent(in)::HminuseS(n*(n+1)/2) ! Hermitian Matrix upper triangle packed
	complex(8),intent(out):: da(n) ! vectors size n
	complex(8)::invHmineS(n)
	complex(8):: norm
	real(8)::rnorm
	integer:: i 
	complex(8) ::zdotc
    external zdotc
	call getinvdiagonalofpacked(n,HminuseS,invHmineS)
write(444,*)"invHmineS	",invHmineS	
write(666,*)"r",residual
	do i=1,n
		da(i)=-invHmineS(i)*residual(i)
	end do
	norm=(0,0)		
    norm=zdotc(n,da,1,da,1)
#ifdef DEBUG
	write(771,*)"da" ,da
#endif
	rnorm=sqrt(dble(norm))
	call zscal(n, DCMPLX(1.0/rnorm,0),da,1)
#ifdef DEBUG

	write(771,*)"da" ,da
	write(772,*)"invHmineS",invHmineS
	write(773)"rnorm",rnorm
	write(774,*)"da" ,da
#endif
end subroutine

subroutine getinvdiagonalofpacked(n,PackedM,Diagonal)
integer,intent(in)::n
complex(8),intent(in)::PackedM(n*(n+1)/2)
complex(8),intent(out)::Diagonal(n)

integer:: diagonalindex,i
diagonalindex=0
do i=1,n
diagonalindex=diagonalindex+i
Diagonal(i)=1.0/PackedM(diagonalindex)
!must be real!
end do

end subroutine
