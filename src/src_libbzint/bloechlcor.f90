!BOP
!
! !ROUTINE: bloechlcor
  
! !INTERFACE:
       subroutine bloechlcor(e,ef,v,blcor)
!      
! !DESCRIPTION:
! 
!   This subroutine calculates the correction to the integration weights for
!   the improved tetrahedron method according to
!   Bl\"ochl {\it et al}, Phys. Rev. B, {\bf 49}, 16223 (1994). The energies at the 
!   corner of the tetrahedron have to be given in increasing order.

! !USES:
             
       implicit none

! !INPUT PARAMETERS:
 
       real(8), intent(in) :: e(4) ! Eigenenergies at the corners of the
!                                     tetrahedron.
 
       real(8), intent(in) :: ef   ! the fermi energy 
       
       real(8), intent(in) :: v    ! Volume of the tetrahedron.

! !OUTPUT PARAMETERS:
 
       real(8), intent(out) :: blcor(4) ! The correction to the 
!                                         integration weights      
       
!   
! !REVISION HISTORY:
!
!   Created: 4th. March 2004, by RGA

! !LOCAL VARIABLES:
 
       integer(4) :: i,j,jm4,index
       real(8)    :: dte
       real(8)    :: sde
       real(8), external :: dos1t
 
 ! !SYSTEM ROUTINES:
 
       intrinsic mod
!EOP
! 
!BOC 
 
      dte=dos1t(e,ef,v)/4.0d+1
      blcor(1:4)=0.0d0
      if(e(4).le.ef) then
        index=4
      else
        index=0
        do i=1,4
          if(e(i).le.ef)index=index+1
        enddo  
      endif
      select case(index)
      case(0,4)
        continue
      case(1,2,3)  
        do i=1,4
         sde=0.0d0
         do j=i,i+2
          jm4=mod(j,4)+1
          sde=sde +e(jm4)
         enddo
         sde=sde-3.0d0*e(i)
         blcor(i)=sde*dte
        enddo
      end select
  
      end subroutine bloechlcor
!EOC
 


