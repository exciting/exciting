
      module gtil
      
      real(8), allocatable :: tilg(:)        ! The tildeg coefficients

      end module gtil
            
!BOP
!
! !ROUTINE: calctildeg
!
! !INTERFACE:
      subroutine calctildeg(lmax)

! !DESCRIPTION:
!
!Calculates $\tilde{g}_{lm,l'm'}$ according to:
!
!\begin{equation}\label{calctilg}
!\tilde{g}_{lm,l'm'}=\sqrt{4\pi}(-1)^{l}\sqrt{\frac{(l+l'+m+m')!(l+l'-m-m')!}%
!{(2l+1)(2l'+1)[2(l+l')+1]!(l+m)!(l-m)!(l'+m')!(l'-m')!}}
!\end{equation}
!
!needed for the calculation of the structure constants, for l and l' = 0
!... lmax and stores them in memory.
!
! !USES:

      use gtil
       
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: lmax
       
! !LOCAL VARIABLES:
       
      integer(4) :: l1,l2,m1,m2,tsize,i

 
! !EXTERNAL ROUTINES: 
      
      real(8), external :: tildeg

!
! !REVISION HISTORY:
!
! Created May 2004 by RGA
! Last modified 30. June 2004 by RGA
!
!EOP
!BOC    
       
      tsize=(lmax+1)*(lmax+2)*(lmax+3)*(3*lmax+2)/12
      if (allocated(tilg)) deallocate(tilg)
      allocate(tilg(tsize))
       
      i=0 
      do l1=0,lmax
        do l2=0,l1
          do m1=-l1,l1
            do m2=0,l2
              i=i+1
              tilg(i)=tildeg(l1,l2,m1,m2)
!              write(96,*)l1,l2,m1,m2,tilg(i)
            enddo ! m2 
          enddo ! m1
        enddo ! l2
      enddo ! l1
      
      end subroutine calctildeg
!EOC
!
!BOP
!
! !ROUTINE: gettildeg
!
! !INTERFACE:
      real(8) function gettildeg(l1,l2,m1,m2)
      
! !DESCRIPTION:

! Gets the value of $\tilde{g}_{l1,m1,l2,m2}$ from the vector
!\verb"tilg"(\verb"i"), where:
!
!\begin{equation}
!\tilde{g}_{l_1m_1,l_2m_2}=(-1)^{l_1-l_2}\texttt{tilg(i)}
!\end{equation}
!
!being
!\begin{equation}
!\texttt{i}=\frac{1}{12}j_1(j_1+1)(j_1+2)(3j_1-1)+j_2(j_2+1)j_1+\frac{1}{2}j_2(j_2-1)%
!+(j_2+1)(m'_1+j_1)+m'_2+j_2+1
!\end{equation}
!
!where we defined:
!
!\begin{equation}
!\begin{aligned}
! j_1=&\textrm{max}(l_1,l_2)\\
! j_2=&\textrm{min}(l_1,l_2)\\
! m'_1=&\textrm{sgn}(m_2) m_1\\
! m'_2=&|m_2|
!\end{aligned}
!\end{equation}
!
! !USES:
 
      use gtil
      
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: l1
      integer(4), intent(in) :: l2
      integer(4), intent(in) :: m1
      integer(4), intent(in) :: m2
      
! !LOCAL VARIABLES:

      integer(4) :: j1,j2,mj1,mj2
      integer(4) :: index1, index2, index3, index4
      integer(4) :: par,tind
      real(8) :: fact


! !INTRINSIC ROUTINES: 


      intrinsic mod
      
! !REVISION HISTORY:
! 
! Creeated May 2004 by RGA
! Last modified: 30. June 2004 by RGA
!
!EOP
!BOC     
      par=mod(abs(l1-l2),2)
!
!     set the indexes j1, j2, mj1, mj2, and the multiplication factor
!
      if(l1.lt.l2)then
        j1=l2
        mj1=m2
        j2=l1
        mj2=m1
!        
!       set the multiplication factor to (-1)^(l1-l2)    
! 
        fact=-2.0d0*par+1.0d0
      else
        j1=l1
        mj1=m1
        j2=l2
        mj2=m2
        fact=1.0d0
      endif
      
      if(mj2.lt.0)then
        mj1=-mj1
        mj2=-mj2
      endif
!
!     get the position of the value searched in the vector tilg
!  
      index1=mj2+j2+1
      index2=(j2+1)*(mj1+j1)
      index3=j2*(j2+1)*j1+j2*(j2-1)/2
      index4=j1*(j1+1)*(j1+2)*(3*j1-1)/12
      tind=index1+index2+index3+index4
      
      gettildeg=fact*tilg(tind)
      
      end function gettildeg

!EOC
