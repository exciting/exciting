      subroutine setpatrd(npar,z,f,a)
!
! this subroutine calculate coeffients in Pade's approximation 
! using Thiele's reciprocal difference method as described in 
!   H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977) 
! or
!   K.-H. Lee and K. J. Chang, Phys. Rev. B 54 R8285 (1996) 
!
! Input:
!  npar       --- order of the Pade's approximation, the number of poles is equal to npar/2
!  z(1..npar) --- complex points, on which, the values of the function to be fitted are given
!  f(1..npar)
!  Output: 
!  a(1..npar) --- coefficients of Pade's approximation
! 
      implicit none 
      integer(4),intent(in)::npar
      complex(8),intent(in)::z(npar),f(npar)
      complex(8),intent(out)::a(npar) 

      integer(4)::n,p
      complex(8)::g(npar,npar)

      g=0.d0
      g(1:npar,1)=f(1:npar) 

      do p=2,npar
        do n=p,npar
          g(n,p)=(g(p-1,p-1)-g(n,p-1))/((z(n)-z(p-1))*g(n,p-1))
        enddo 
      enddo 
      do n=1,npar
        a(n)=g(n,n)
      enddo 
      end subroutine 
