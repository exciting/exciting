!BOP
!
! !ROUTINE: stweight_itaylor
!
! !INTERFACE:

      subroutine stweight_itaylor(deltae_vert,omeg,weight_vert)
!
! !DESCRIPTION:
!
! This subroutine calculates the weight on the whole small tetrahedron
! in which the $k$ states are fully occupied and $k-q$ states are fully 
! unoccupied. This is for the $sigfreq=3$ case when we consider the 
! imaginary frequency. 
!                                             
! !INPUT PARAMETERS:
      implicit none
      
      real(8), intent(in) :: deltae_vert(4)        ! difference of the energy
                                                   ! in k-mesh tetrahedron vertices 
                                                   ! and k-q mesh tetrahedron vertices.
      real(8), intent(in) :: omeg                  ! the frequency omega to be calculated
 
! !OUTPUT PARAMETERS:            

      real(8), intent(out) :: weight_vert(4)! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: ivert,j,k

      real(8)    :: omt2,omt4,denom1,denom3,w01,n03,w03

      real(8), dimension(4) :: ev


 
! !REVISION HISTORY:
!
! Created 01.02.2006 by RGA.
!

!EOP
!BOC
      omt2=omeg*omeg
      omt4=omt2*omt2
      denom1=6.0d+1*omt2
      denom3=4.2d+2*omt4
      do ivert=1,4
        do j=1,4
          k=mod(j+ivert-2,4)+1
          ev(j)=deltae_vert(k)
        enddo
        w01=-(2.0d0*ev(1)+ev(2)+ev(3)+ev(4))/denom1
        n03=4.0d0*ev(1)**3+3.0d0*ev(1)**2*(ev(2)+ev(3)+ev(4))+          &
     &      2.0d0* ev(1)*(ev(2)**2+ev(3)**2+ev(4)**2+ev(2)*ev(3)+ev(2)* &
     &      ev(4)+ev(3)*ev(4))+ev(2)**3+ev(3)**3+ev(4)**3+ev(2)**2*     &
     &      ev(3)+ev(2)**2*ev(4)+ev(3)**2*ev(4)+ev(2)*ev(3)**2+ev(2)*   &
     &      ev(4)**2+ev(3)*ev(4)**2+ev(2)*ev(3)*ev(4)
        w03=n03/denom3
        weight_vert(ivert)=w01+w03
      enddo
      return
      
      end subroutine stweight_itaylor
!EOC      
