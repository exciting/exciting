
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

      subroutine edifwtaylor(v,omeg,wt)
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
      do i=1,4
        do j=1,4
          ind=mod(j+i-2,4)+1
          ev(j)=v(ind)
        enddo
        w01=-(2.0d0*ev(1)+ev(2)+ev(3)+ev(4))/denom1
        n03=4.0d0*ev(1)**3+3.0d0*ev(1)**2*(ev(2)+ev(3)+ev(4))+           &
     &      2.0d0* ev(1)*(ev(2)**2+ev(3)**2+ev(4)**2+ev(2)*ev(3)+ev(2)* &
     &      ev(4)+ev(3)*ev(4))+ev(2)**3+ev(3)**3+ev(4)**3+ev(2)**2*     &
     &      ev(3)+ev(2)**2*ev(4)+ev(3)**2*ev(4)+ev(2)*ev(3)**2+ev(2)*   &
     &      ev(4)**2+ev(3)*ev(4)**2+ev(2)*ev(3)*ev(4)
        w03=n03/denom3
        wt(i)=w01+w03
      enddo
      return
      
      end subroutine edifwtaylor
!EOC      
