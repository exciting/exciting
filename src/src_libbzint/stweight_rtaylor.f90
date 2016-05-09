!BOP
!
! !ROUTINE: stweight_rtaylor
!
! !INTERFACE:

      subroutine stweight_rtaylor(deltae_vert,freq,weight_vert)
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
!                                    in k-mesh tetrahedron vertices 
!                                  and k-q mesh tetrahedron vertices.

      real(8), intent(in) :: freq ! the frequency omega to be calculated

 
! !OUTPUT PARAMETERS:            

      real(8), intent(out) :: weight_vert(4)! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: ivert,j,k,isign_om

      real(8) :: denom0, denom1, denom2, denom3, denom4
   
      real(8) :: w00, w01, w02, w03
      real(8) :: omeg

      real(8), dimension(4) :: ev
      real(8), dimension(2,4) :: weight_tmp


 
! !REVISION HISTORY:
!
! Created 01.02.2006 by XZL.
!

!EOP
!BOC
      do isign_om=1,2
        omeg=dble(3-2*isign_om)*freq
        denom0=omeg
        denom1=omeg**2
        denom2=omeg**3
        denom3=omeg**4
        denom4=omeg**5
        do ivert=1,4
          do j=1,4
            k=mod(j+ivert-2,4)+1
            ev(j)=deltae_vert(k)
          enddo
          w00=1.0d0/(2.4d+1*denom0)
          w01=(ev(1)+2.0d0*ev(2)+ev(3)+ev(4))/(1.2d+2*denom1)
          w02=(ev(1)**2+2.0d0*ev(1)*ev(2)+3.0d0*ev(2)**2+ev(1)*ev(3)+     &
     &        2.0d0*ev(2)*ev(3)+ev(3)**2+ev(1)*ev(4)+2.0d0*ev(2)*ev(4)+   &
     &        ev(3)*ev(4)+ev(4)**2)/(3.6d+2*denom2)
          w03=(ev(1)**3+2.0d0*ev(1)**2*ev(2)+3.0d0*ev(1)*ev(2)**2+        &
     &        4.0d0*ev(2)**3+ev(1)**2*ev(3)+2.0d0*ev(1)*ev(2)*ev(3)+      &
     &        3.0d0*ev(2)**2*ev(3)+ev(1)*ev(3)**2+2.0d0*ev(2)*ev(3)**2+   &
     &        ev(3)**3+ev(1)**2*ev(4)+2.0d0*ev(1)*ev(2)*ev(4)+            &
     &        3.0d0*ev(2)**2*ev(4)+ev(1)*ev(3)*ev(4)+2.0d0*ev(2)*ev(3)*   &
     &        ev(4)+ev(3)**2*ev(4)+ev(1)*ev(4)**2+2.0d0*ev(2)*ev(4)**2+   &
     &        ev(3)*ev(4)**2+ev(4)**3)/(8.4d+2*denom3)
          weight_tmp(isign_om,ivert)=w00+w01+w02+w03
        enddo ! ivert
      enddo ! isign_om
      do ivert=1,4
        weight_vert(ivert)=weight_tmp(1,ivert)+weight_tmp(2,ivert)
      enddo  

      return
   
      end subroutine stweight_rtaylor
!EOC      
