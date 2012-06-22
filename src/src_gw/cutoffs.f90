!BOP
!
! !MODULE: incgamma
      module incgamma
      
! !LOCAL VARIABLES:
   
      real(8), private :: etx ! exp(-x), declared here because it is used
!                               by all the calls of the function, to avoid
!                               recalculating it every time.

      real(8), private :: sqx ! sqrt(x), declared here because it is used
!                               by all the calls of the function, to avoid
!                               recalculating it every time.

! !INTERFACE:

      interface gammaincc
        module procedure gammaincc_int
        module procedure gammaincc_hin
      end interface gammaincc

!
! !DESCRIPTION:
!
! This module calculates the Incomplete gamma function, defined by:
!
!\begin{equation}
!\Gamma(a,x)=\int\limits_x^{\infty}e^{-t} t^{a-1}dt
!\end{equation}
! \noindent 
! for values of a integer and halfinteger recursively,
! according to the equation:
!
!\begin{equation}
! \Gamma(a+1,x)= a\Gamma(a,x)+x^ae^{-x}
!\end{equation}
!
!!USAGE: 
!
!      result=incgam(n,x)
!
!  where:
!
!\begin{itemize}
! \item n (real(8)): order of the gamma function, the number must be
! integer or half integer
! \item x (real(8)): point where the gamma function is to be
! calculated. x must be positive
!\end{itemize}
!
! !EXTERNAL ROUTINES: 
!
!\begin{itemize}
! \item \verb"e1xb": calculates the exponential integral function
! \item \verb"erf": calculates the complementary error function. Depending
! on the compiler this can also be an internal procedure (ifc 6.0
! has  it). If not, the file erf.f90 is included in the package.
!\end{itemize}
!
! !PUBLIC MEMBER FUNCTIONS:
!
! real(8) :: incgam(n,x)
!
!EOP
      contains
!BOP
!
! !IROUTINE: gammaincc_int
! 
! !INTERFACE:
      recursive function gammaincc_int(n,x) result(gmi)
      
! !INPUT PARAMETERS:      
        implicit none
        
        real(8) :: x ! value at which the gamma function is 
!                      calculated
        
        integer(4) :: n ! order of the gamma function

! !OUTPUT PARAMETERS:        
        
        real(8) :: gmi ! value of the gamma function
!        
! !DESCRIPTION:
!
! Calculates the incomplete gamma function of integer order
!
!EOP
!
!BOC
        real(8), external :: e1xb
        if(n.eq.0)then
          gmi=e1xb(x)
        elseif(n.eq.1)then
          gmi = etx
        else
          gmi = (x**(n-1))*etx+dble(n-1)*gammaincc(n-1,x)
        endif
      end function gammaincc_int
!EOC
!
!BOP
!
! !IROUTINE: gammaincc_hin
! 
! !INTERFACE:
      recursive function gammaincc_hin(in,x,n) result(gmh)
! !INPUT PARAMETERS:      
        implicit none
        
        integer(4) :: in ! order of the gamma function - 1/2

        real(8) :: x ! value at which the gamma function is calculated
        
        real(8) :: n ! order of the gamma function

! !OUTPUT PARAMETERS:        
        
        real(8) :: gmh ! value of the gamma function
!        
! !DESCRIPTION:
!
! Calculates the incomplete gamma function of halfinteger order
!
!EOP
!
!BOC
        real(8), external :: derfc
        real(8) :: pi
        if(in.eq.0)then
          pi = 4.0d0*datan(1.0d0)
          gmh=dsqrt(pi)*derfc(sqx)
        else
          gmh = (x**(in-1))*sqx*etx+(n-1.0d0)*&
     &         gammaincc(in-1,x,n-1.0d0)
        endif
      end function gammaincc_hin
!EOC
!
!BOP
! 
! !IROUTINE: incgam
!
! !INTERFACE:
      function incgam(n,x)
      
! !INPUT PARAMETERS:
        implicit none
        
        real(8) :: x     ! Argument of the gamma function

        real(8) :: n     ! order of the gamma function

! !OUTPUT PARAMETERS:

        real(8) :: incgam ! value of the gamma function (result)

! !DESCRIPTION:
!
! The function works as an interface, it checks whether the index
! is integer or half integer and calls the corresponding procedure
! within the module
!
! !LOCAL VARIABLES:

        integer(4) :: in ! integer part of n
        
! !REVISION HISTORY:
!
! Created Jan. 2004        
!        
!EOP
!
!BOC
        etx = dexp(-1.0d+0*x)
        in=idint(n)
        if(n-dble(in).le.0.25)then
          incgam=gammaincc(in,x)
        else
          sqx=dsqrt(x)
          incgam=gammaincc(in,x,n)
        endif
      end function incgam
!EOC
      end module incgamma

!BOP
!
! !ROUTINE: rcutoff
!
! !INTERFACE:
      function rcutoff(tol,eta,lambdamax) result(rcf)

! !DESCRIPTION:
!      
!  Estimates the cutoff radius of the sums in real space for the
! calculation of the structure constants by the solving the equation:
!
!\begin{equation}
!\mathfrak{E}_{R,\lambda}^{\textrm{tol}}=\left\{%
!\begin{array}{ll}
!\frac{4\pi}{(\lambda -2)\Gamma(\lambda+\tfrac{1}{2})}%
!\left(\frac{\Gamma[\tfrac{\lambda}{2}+\tfrac{3}{2},\left(\tfrac{R_c}{\eta}\right)^2]}%
!{\eta^{\lambda-2}}-\frac{\Gamma[\lambda+\tfrac{1}{2},\left(\tfrac{R_c}{\eta}\right)^2]}%
!{R_c^{\lambda-2}}\right)&\lambda \neq 2\\
!\frac{4\pi}{\Gamma(\tfrac{5}{2})}\left[\tfrac{\eta}{R_c}%
!\Gamma[3,\left(\tfrac{R_c}{\eta}\right)^2]-\Gamma[\tfrac{5}{2},\left(\tfrac{R_c}{\eta}\right)^2]\right]&
!\lambda=2\\
!\end{array}
!\right.
!\end{equation}
!
! and taking the maximum value of $R_c$ obtained for $\lambda = 1...$
!\verb"lambdamax".      
!
! !USES:
      
      use modmain,  only: pi
      use modgw,    only: avec
      use incgamma, only: incgam
      
! !INPUT PARAMETERS:

      implicit none

      real(8),    intent(in) :: tol ! The tolerance for the 
!                                     convergence of the lattice sum

      real(8),    intent(in) :: eta 
 
      integer(4), intent(in) :: lambdamax

! !OUTPUT PARAMETERS:      
 
      real(8) :: rcf ! The maximum cutoff radius

! !LOCAL VARIABLES:

      integer(4) :: l1
      integer(4) :: i

      real(8) :: rl
      real(8) :: x
      real(8) :: gmm
      real(8) :: aleng(3)
      real(8), allocatable :: eps(:)
      real(8) :: rnot
      real(8) :: gaml12
      real(8) :: gaml32
      real(8) :: prefac
      real(8), allocatable :: rct(:)

      real(8), external :: higam

!
! !REVISION HISTORY:
! 
! Created: January 2004 by RGA
! Last Modified: 6. August 2004 by RGA
! Revisited: 23 May 2011 by DIN
!
!EOP
!BOC
      allocate(rct(lambdamax+1))
      allocate(eps(lambdamax+1))
      do i=1,3
        aleng(i)=sqrt(avec(1,i)*avec(1,i)+ &
       &              avec(2,i)*avec(2,i)+ &
       &              avec(3,i)*avec(3,i))
      enddo
      rnot=maxval(aleng)
      rct = 5.0d+1
      do i=1,100
        x = 5.0d-1*dble(i)
        do l1=0,lambdamax
          if(l1.ne.2)then
            rl =5.0d-1*(dble(l1+1))
            gaml32=incgam(rl,x*x)
            rl = dble(l1)+5.0d-1
            gaml12=incgam(rl,x*x)
            gmm = higam(l1)
            prefac = 4.0d0*pi/(dble(l1-2)*gmm)
            eps(l1+1)=dabs(prefac*(gaml32-gaml12/(x**(l1-2)))/        &
     &                (eta**(l1-2)))
            if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
          else
            gaml32=incgam(3.0d0,x*x)
            gaml12=incgam(2.5d0,x*x)
            gmm = higam(2)
            prefac = 4.0d0*pi/gmm
            eps(l1+1)= dabs(prefac*(gaml32/x-gaml12))
            if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
          endif
        enddo ! l1
!        write(52,10)x,(eps(l1),l1=1,lambdamax+1)
      enddo ! i
!      do l1=1,lambdamax+1
!        write(52,*)'RGA: rcutoff:',l1-1,rct(l1)
!      enddo
      rcf=maxval(rct)*eta
!      write(52,*)'max. real cutoff =',rcf
      
      deallocate(rct)
      deallocate(eps)
      
   10 format(12(e14.7e3,1x))
   
      end function rcutoff
!EOC

!BOP
!
! !ROUTINE: gcutoff
!
! !INTERFACE:
      function gcutoff(tol,eta,lambdamax) result(rcf)
      
! !DESCRIPTION:
!      
!  Estimates the cutoff radius of the sums in reciprocal space for the
! calculation of the structure constants by the solving the equation:
!
!\begin{equation}
!\mathfrak{E}_{G,\lambda}^{\textrm{tol}}=\frac{8(\pi)^{\frac{5}{2}}}{\Omega%
!\Gamma(\lambda+\tfrac{1}{2})\eta^{\lambda+1}}%
!\Gamma\left[\tfrac{\lambda+1}{2},\left(\tfrac{\eta G_c}{2}\right)^2\right]
!\end{equation}
!
! and taking the maximum value of $G_c$ obtained for $\lambda = 1...$
!\verb"lambdamax".      
!
! !USES:
      
      use modmain,   only: pi,omega,bvec
      use incgamma,  only: incgam

! !INPUT PARAMETERS:

      implicit none

      real(8),    intent(in) :: tol ! The tolerance for the 
!                                     convergence of the lattice sum
      real(8),    intent(in) :: eta 
      integer(4), intent(in) :: lambdamax

! !OUTPUT PARAMETERS:      
 
      real(8) :: rcf ! The maximum cutoff radius

! !LOCAL VARIABLES:

      integer(4) :: l1
      integer(4) :: i

      real(8) :: gaml32
      real(8) :: gmm
      real(8) :: prefac
      real(8) :: rl
      real(8) :: rnot
      real(8) :: x
      real(8) :: bleng(3)
      real(8), allocatable :: eps(:)
      real(8), allocatable :: rct(:)

 
! !EXTERNAL ROUTINES: 


      real(8), external :: higam
!
! !REVISION HISTORY:
! 
! Created: January 2004 by RGA
! Last Modified: 6. August 2004 by RGA
!
!EOP
!BOC

      allocate(rct(lambdamax+1))
      allocate(eps(lambdamax+1))
      do i=1,3
        bleng(i)=sqrt(bvec(1,i)*bvec(1,i)+ &
       &              bvec(2,i)*bvec(2,i)+ &
       &              bvec(3,i)*bvec(3,i))
      enddo
      rnot=maxval(bleng)
      rct(:) = 5.0d+1
      do i=1,100
        x = 5.0d-1*dble(i)
        do l1=0,lambdamax
          rl =5.0d-1*(dble(l1+1))
          gaml32=incgam(rl,x*x)
          gmm = higam(l1)
          prefac = 8.0d0*pi*pi*dsqrt(pi)/(omega*gmm*eta**(l1+1))
          eps(l1+1)=dabs(prefac*gaml32)
          if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
        enddo ! l1
!        write(52,10)x,(eps(l1),l1=1,lambdamax+1)
      enddo ! i
!      do l1=1,lambdamax+1
!        write(52,*)'RGA: gcutoff:',l1-1,rct(l1)
!      enddo
      rcf=maxval(rct)*2.0d0/eta
!      write(52,*)'max. reciprocal cutoff =',rcf
      
   10 format(12(e14.7e3,1x))
   
      deallocate(rct)
      deallocate(eps)

      end function gcutoff
!EOC      
