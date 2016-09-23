      subroutine acpatrd(npar,z,zn,apar,fz,dfz) 
!
! this subroutine calculate fitting value by Pade's approximation
! using Thiele's reciprocal difference method as described in
!   H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977)
! or
!   K.-H. Lee and K. J. Chang, Phys. Rev. B 54 R8285 (1996)
!
! Input:
!  rez, imz   --- real and imaginary part of z
!  npar       --- order of the Pade's approximation, the number of poles is equal to npar/2
!  zn(1..npar) --- complex points, on which, the values of the function to be fitted are given
!  a(1..npar)  --- PA parameters calculated in the subroutine inipatrd
!  Output:
!   fz       --- f(z) 
!   dfz      --- f'(z)
!


      implicit none 
      integer(4),intent(in)::npar
      complex(8),intent(in)::z,zn(npar),apar(npar) 
      complex(8),intent(out)::fz,dfz

      integer(4)::n
      complex(8)::coef
      complex(8)::aa(0:npar),bb(0:npar),daa(0:npar),dbb(0:npar) 

      
      aa(0)=cmplx(0.d0,0.d0)
      aa(1)=apar(1)
      bb(0)=cmplx(1.d0,0.d0)
      bb(1)=bb(0)
      daa(0:1)=cmplx(0.d0,0.d0)
      dbb(0:1)=cmplx(0.d0,0.d0)

      do n=2,npar
        coef=(z-zn(n-1))*apar(n)
        aa(n)=aa(n-1)+coef*aa(n-2)
        bb(n)=bb(n-1)+coef*bb(n-2)
        daa(n)=daa(n-1)+coef*daa(n-2)+apar(n)*aa(n-2)
        dbb(n)=dbb(n-1)+coef*dbb(n-2)+apar(n)*bb(n-2)
      enddo 

      fz=aa(npar)/bb(npar)
      dfz=(daa(npar)*bb(npar)-aa(npar)*dbb(npar))/(bb(npar)*bb(npar))
      end subroutine 
