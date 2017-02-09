!BOP
!
! !ROUTINE: average_degen_convw 
!
! !INTERFACE:
       subroutine average_degen_convw(cweight)

!     
! !DESCRIPTION:
!
!   This subroutine averages the weights of degenerate states, an
! intrinsic error of the tetrahedron method.
!  

! !USES:
 
       use tetra_internal
       
       implicit none      
       
! !INPUT/OUTPUT PARAMETERS:
 
       real(8), intent(inout)   :: cweight(nband,nband,nirkp) ! the value of the integral
       
! !LOCAL VARIABLES:

      integer(4) :: ik
      integer(4) :: jk
      integer(4) :: ib
      integer(4) :: ib1
      integer(4) :: jb
      integer(4) :: degeneracy

      real(8) :: edist
      real(8) :: avg_weight  
      real(8) :: sum_weight
      real(8), parameter :: zerotol = 1.0d-3   
      
      logical :: are_deg        
       
! !INTRINSIC ROUTINES:
       
       
! !EXTERNAL ROUTINES:

       
! !REVISION HISTORY:
!
!   Created: 4th. March 2004 by RGA
!
!EOP
!BOC

      do ik = 1, nirkp
        ib = 1
        do while (ib.lt.nband)
          degeneracy = 0
          are_deg = .true.
          do while (are_deg)
            degeneracy = degeneracy +1
            ib1 = ib + degeneracy
            if(ib1.le.nband)then ! avoids running out argument
! Test if the states are degenerate, then are_deg = true              
              edist=abs(eband(ib1,ik)-eband(ib,ik))
              are_deg = edist.le.zerotol
            else
              are_deg = .false.
            endif    
          enddo

          do jb = 1, nband 
!
! Average the weights of the degenerate states
!          
            sum_weight = 0.0d0
            do ib1 = ib, ib + degeneracy - 1
              sum_weight = sum_weight + cweight(ib1,jb,ik)
            enddo
            avg_weight = sum_weight / dble(degeneracy)
!
! Reasign the average weight to all the degenerate states
!
            do ib1 = ib, ib + degeneracy - 1
              cweight(ib1,jb,ik) = avg_weight
            enddo
          enddo  
!
! Jump to the next non-degenerate states
!
            ib=ib+degeneracy
        enddo
      enddo    
!-----------------------------------------------------
!             Now the same for the second index
!-----------------------------------------------------
      do ik = 1, nirkp
        jk=klinkq(ik)
        ib = 1
        do while (ib.lt.nband)
          degeneracy = 0
          are_deg = .true.
          do while (are_deg)
            degeneracy = degeneracy +1
            ib1 = ib + degeneracy
            if(ib1.le.nband)then ! avoids running out argument
! Test if the states are degenerate, then are_deg = true              
              edist=abs(eband(ib1,jk)-eband(ib,jk))
              are_deg = edist.le.zerotol
            else
              are_deg = .false.
            endif    
          enddo

          do jb = 1, nband 
!
! Average the weights of the degenerate states
!          
            sum_weight = 0.0d0
            do ib1 = ib, ib + degeneracy - 1
              sum_weight = sum_weight + cweight(jb,ib1,ik)
            enddo
            avg_weight = sum_weight / dble(degeneracy)
!
! Reasign the average weight to all the degenerate states
!
            do ib1 = ib, ib + degeneracy - 1
              cweight(jb,ib1,ik) = avg_weight
            enddo
          enddo  
!
! Jump to the next non-degenerate states
!
            ib=ib+degeneracy
        enddo
      enddo    

      end subroutine average_degen_convw
!EOC      
            
           
