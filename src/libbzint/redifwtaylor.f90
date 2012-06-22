
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: edifwtaylor
!
! !INTERFACE:

      subroutine redifwtaylor(v,omeg,wt)
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
      
      real(8), intent(in) :: v(4)        ! difference of the energy
!                                    in k-mesh tetrahedron vertices 
!                                  and k-q mesh tetrahedron vertices.

      real(8), intent(in) :: omeg ! the frequency omega to be calculated

 
! !OUTPUT PARAMETERS:            

      real(8), intent(out) :: wt(4)! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: i,j,ind

      real(8) :: denom0, denom1, denom2, denom3, denom4
   
      real(8) :: w00, w01, w02, w03

      real(8), dimension(4) :: ev


 
! !REVISION HISTORY:
!
! Created 01.02.2006 by XZL.
!

!EOP
!BOC
      denom0=omeg
      denom1=omeg**2
      denom2=omeg**3
      denom3=omeg**4
      denom4=omeg**5
      do i=1,4
        do j=1,4
          ind=mod(j+i-2,4)+1
          ev(j)=v(ind)
        enddo
        w00=1.0d0/(2.4d+1*denom0)
        w01=(ev(1)+2.0d0*ev(2)+ev(3)+ev(4))/(1.2d+2*denom1)
        w02=(ev(1)**2+2.0d0*ev(1)*ev(2)+3.0d0*ev(2)**2+ev(1)*ev(3)+     &
     &      2.0d0*ev(2)*ev(3)+ev(3)**2+ev(1)*ev(4)+2.0d0*ev(2)*ev(4)+   &
     &      ev(3)*ev(4)+ev(4)**2)/(3.6d+2*denom2)
        w03=(ev(1)**3+2.0d0*ev(1)**2*ev(2)+3.0d0*ev(1)*ev(2)**2+        &
     &      4.0d0*ev(2)**3+ev(1)**2*ev(3)+2.0d0*ev(1)*ev(2)*ev(3)+      &
     &      3.0d0*ev(2)**2*ev(3)+ev(1)*ev(3)**2+2.0d0*ev(2)*ev(3)**2+   &
     &      ev(3)**3+ev(1)**2*ev(4)+2.0d0*ev(1)*ev(2)*ev(4)+            &
     &      3.0d0*ev(2)**2*ev(4)+ev(1)*ev(3)*ev(4)+2.0d0*ev(2)*ev(3)*   &
     &      ev(4)+ev(3)**2*ev(4)+ev(1)*ev(4)**2+2.0d0*ev(2)*ev(4)**2+   &
     &      ev(3)*ev(4)**2+ev(4)**3)/(8.4d+2*denom3)
        ind=mod(i,4)+1
        wt(ind)=w00+w01+w02+w03
      enddo
      return
      
   
      end subroutine redifwtaylor
!EOC      
