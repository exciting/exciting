      subroutine getsac(iop,nomeg,npar,en,ein,omega,apar,sc,dsc) 
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
      integer(4),intent(in)::iop,nomeg,npar
      real(8),intent(in)::en
      real(8),intent(in)::omega(nomeg)
      complex(8),intent(in)::ein,apar(npar) 
      complex(8),intent(out)::sc,dsc

      complex(8)::comega(npar)
      integer(4)::ip,iw,iwpa(npar)

      
      

      if(iop.eq.1) then 
        call acrgn(npar,ein,apar,sc,dsc)
      elseif(iop.eq.0) then   
        call setwpa(nomeg,npar,omega,iwpa)
        do ip=1,npar
          iw=iwpa(ip)
          comega(ip)=cmplx(0.d0,dsign(omega(iw),en))
        enddo 
        call acpatrd(npar,ein,comega,apar,sc,dsc)
      elseif(iop.eq.2) then  
         call setwpa(nomeg,npar,omega,iwpa)
         do ip=1,npar
           iw=iwpa(ip)
           comega(ip)=cmplx(-en,dsign(omega(iw),en))
         enddo
         call acpatrd(npar,ein-en,comega,apar,sc,dsc)
      endif 




      end subroutine 
