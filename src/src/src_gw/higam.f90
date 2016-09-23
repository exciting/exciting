!BOP
! !ROUTINE: higam
! 
! !INTERFACE: 
        recursive function higam(n) result(tgam)
! !DESCRIPTION:
!
! This function calculates the gamma function for half integer values
! $\Gamma(n+\tfrac{1}{2})$ recursively, using the relation:
! $\Gamma(n+\tfrac{1}{2})=\frac{2n-1}{2}\Gamma(n-\tfrac{1}{2})$, and
! $\Gamma(\tfrac{1}{2})=\sqrt{\pi}$.
!
! !INPUT PARAMETERS:        
        implicit none
        integer(4) :: n
        
! !OUTPUT PARAMETERS:        
        real(8) :: tgam
        
! !DEFINED PARAMETERS:        
        real(8) :: pi
        pi = 4.0d+0*datan(1.0d+0)
!EOP
!BOC        
        if(n.lt.0)then
          write(*,*)"rga: higam: n mus be positive"
          stop
        elseif(n.eq.0)then
          tgam = dsqrt(pi)
        else
          tgam= 5.0d-1*dble(2*n-1)*higam(n-1)
        endif
        end function higam
!EOC
