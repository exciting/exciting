subroutine calcupdatevector(n,residual,HminuseS,da)
	integer, intent(in)::n
    complex(8),intent(in)::residual(n)
	complex(8),intent(in)::HminuseS(n*(n+1)/2) ! Hermitian Matrix upper triangle packed
	complex(8),intent(out):: da(n) ! vectors size n
	complex(8)::invHmineS(n)
	real(8):: norm
	integer:: i
	call getinvdiagonalofpacked(n,HminuseS,invHmineS)
write(444,*)"invHmineS	",invHmineS	
write(666,*)"r",residual
	do i=1,n
		da(i)=-invHmineS(i)*residual(i)
	end do
	write(777,*)"da" ,da
	norm=0
	do i=1,n
		norm=(norm+dble(conjg(da(i))*da(i)))
	end do
	write(555,*)"da norm",norm
	norm=sqrt(norm)
	write(555,*)"da norm",norm
	call zscal(n, DCMPLX(1.0/norm),da,1)
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
