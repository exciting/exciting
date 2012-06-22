!BOP
!
! !ROUTINE: getdjmm
!
! !INTERFACE:
      complex(8) function getdjmm(isym,l,m1,m2)
      
! !DESCRIPTION:
!
! Extracts the value of $D^j_{m_1m_2}$ from the djmm vector by using:
!
!\begin{equation}
!D^j_{m_1m_2}=\left\{
!\begin{array}{ll}
!\texttt{djmm(i)} & m_2 \ge 0 \\
!(-1)^{m_1-m_2}\texttt{djmm(i)}^* & m_2 < 0
!\end{array}
!\right.
!\end{equation}
!
!where the index is given by
!
!\begin{equation}
!\texttt{i}=\tfrac{1}{6}l(l+1)(4j-1)+(l+1)(\textrm{sgn}(m_2)m_1+l)+|m_2|+1
!\end{equation}

! !USES:

      use modgw, only: djmm
      
! !INPUT PARAMETERS:
      
      implicit none
      
      integer(4), intent(in) :: isym
      
      integer(4), intent(in) :: l        
      
      integer(4), intent(in) :: m1
      
      integer(4), intent(in) :: m2
      
! !LOCAL VARIABLES:

      integer(4) :: i,par
      
      real(8) :: alpha
      

! !INTRINSIC ROUTINES: 

      intrinsic iabs
      intrinsic sign
      
! !REVISION HISTORY:
!
! Created 11. August 2004      
!EOP
!BOC
!
!     first check that the input parameters are physical
!
      if((iabs(m2).le.l).and.(iabs(m1).le.l))then
!
!       calculate the index of the djmm vector
!       
        i=l*(l+1)*(4*l-1)/6+(l+1)*(sign(1,m2)*m1+l)+iabs(m2)+1
        if(m2.ge.0)then      
!
!         For m2>=0 D^j_m1m2=djmm(i)
!        
          getdjmm=djmm(isym,i)
        else
!
!         For m2<0 D^j_m1m2=(-1)^(m1-m2)*djmm^*(i)
!        
          par=mod(abs(m1-m2),2)
          alpha=1.0d0-2.0d0*dble(par)
          getdjmm=alpha*conjg(djmm(isym,i))

        endif
!
!       Successful exit
!        
        return
      else
!
!       Error handling
!      
        write(*,*)'getdjmm: error'
        write(*,*)'unphysical values as input'
        write(*,*)'l =',l,'m1 =',m1,'m2 =',m2
        stop 'error in getdjmm'
      endif
      end function getdjmm  
!EOC               
