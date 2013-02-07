!BOP
!
! !ROUTINE: genauxf
!
! !INTERFACE: 
      subroutine genauxf(iq,beta,f1,f2)

! !DESCRIPTION:
!
! Given the index \texttt{iq} of $\vec{q}$, this
! subroutine generates the auxiliary functions $F_1(\vec{q})$ and $F_2(\vec{q})$
! according to the formulas:
!
! \begin{subequations}\label{genauxf-01}
! \begin{align}
! F_1(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha|\vec{q}+\vec{G}_i|^2}}{|\vec{q}+\vec{G}_i|}}\\
! F_2(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha|\vec{q}+\vec{G}_i|^2}}{|\vec{q}+\vec{G}_i|^2}}
! \end{align}
! \end{subequations}
! 
! where $\beta=\left(\frac{\Omega}{6\pi^2}\right)^{\frac{1}{3}}$
!
! !USES:

      use modmain
      use modgw      

! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: iq    ! Index of the q-vector
!                                       at iq       
      real(8),    intent(in) :: beta  ! Parameter of the function
      
! !OUTPUT PARAMETERS:

      real(8),    intent(out) :: f1   ! Value of the auxiliary function F_1
!                                       at iq       
      real(8),    intent(out) :: f2   ! Value of the auxiliary function F_2
      
! !LOCAL VARIABLES:

      integer(4) :: ipw              ! (Counter) runs over plane waves.
      integer(4) :: ipwin            ! Initial value of ipw
      real(8) :: modgpq              ! Length of q+G squared
      real(8) :: expagpq             ! exp(beta*|q+g|^2
      
      real(8), dimension(3) :: qvec  ! cartesian coords. of the q point
      real(8), dimension(3) :: gvec  ! cartesian coords. of the vector G.
      real(8), dimension(3) :: gpq   ! cartesian coords. of the G+q.
 
! !INTRINSIC ROUTINES: 

      intrinsic exp
      intrinsic sqrt
      

! !EXTERNAL ROUTINES: 


! !REVISION HISTORY:
!
! Created 23.06.05 by RGA.
! Revisited June 2011 by DIN
!
!EOP
!BOC
!
!     Calculate the cartessian coordinates of the q-point
!      
      qvec(1:3)=vqc(1:3,iq)
      ipwin=1
      if(Gamma)ipwin=2 
!
!     Initializations
!
      f1=0.0d0
      f2=0.0d0
!
!     Loop over G-vectors
!      
      do ipw = ipwin, ngbarc(iq)
!
!       Calculate the cartessian coordinates of G
!      
        gvec(1:3)=vgc(1:3,igqigb(ipw,iq))
!
!       Calculate the length of G+q squared 
!
        gpq(1:3) = gvec(1:3) + qvec(1:3)
        modgpq=gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3)
        expagpq=exp(-beta*modgpq)
!
!       Accumulate the terms into f1 and f2
!
        f1=f1+expagpq/sqrt(modgpq)
        f2=f2+expagpq/modgpq

      enddo ! ipw

      return
      end subroutine genauxf
!EOC                
