!BOP
!
! !ROUTINE: rot1tosph
!
! !INTERFACE:
      subroutine rot1tosph(rotcart,rotsph)
      
! !DESCRIPTION:
!
! Transform a rotation matrix from cartesian $(x,y,z)$ to spherical basis
!$(\tfrac{1}{\sqrt{2}}(x-iy),z,-\tfrac{1}{\sqrt{2}}(x+iy))$, that is
!$(Y_{11},Y_{10},Y_{11})$.
!
! !INPUT PARAMETERS:
 
      implicit none
      
      real(8), intent(in) :: rotcart(3,3) ! the rotation matrix in
!                                              cartesian basis 

! !OUTPUT PARAMETERS:
      
      complex(8), intent(out) :: rotsph(3,3) ! the rotation matrix in 
!                                              spherical basis     

! !LOCAL VARIABLES:
      
      integer(4) :: i, j
      
      complex(8), dimension(3,3) :: umat     ! the transformation matrix
      
      complex(8), dimension(3,3) :: tmat1,tmat2     ! temporary matrices
      
      complex(8), dimension(3,3) :: crcart   ! unitary rotation matrix in
!                                              the cartesian basis
      
! !DEFINED PARAMETERS:

      real(8), parameter :: isq2 = 0.7071067811865475727373109293694142252

!                            1/sqrt(2)
 
!      data umat /(isq2,0.0d0),  (0.0d0,0.0d0), (-isq2,0.0d0),         &
!     &           (0.0d0,-isq2), (0.0d0,0.0d0), (0.0d0,-isq2),         &
!     &           (0.0d0,0.0d0), (1.0d0,0.0d0), (0.0d0,0.0d0)/  
     
      intrinsic mod
      intrinsic conjg
     
!
! !REVISION HISTORY:
!
! Created 10th. August 2004 by (RGA)
!
!EOP
!BOC

!
!        Initialize umat (DIN: gfortran seems to do not support DATA statement)
!
        umat(:,1) = (/ cmplx(isq2,0.0d0),  cmplx(0.0d0,0.0d0), cmplx(-isq2,0.0d0) /)
        umat(:,2) = (/ cmplx(0.0d0,-isq2), cmplx(0.0d0,0.0d0), cmplx(0.0d0,-isq2) /)
        umat(:,3) = (/ cmplx(0.0d0,0.0d0), cmplx(1.0d0,0.0d0), cmplx(0.0d0,0.0d0) /)

!
!       calculate the inverse of umat
!      
        do i=1,3
          do j=1,3
            tmat2(i,j)=conjg(umat(j,i))
          enddo
        enddo    
!
!       calculate the unitary rotation matrix in cartesian basis
!        
        do i=1,3
          do j=1,3
            crcart(i,j)=cmplx(rotcart(i,j),0.0d0,8)
          enddo
        enddo    
        tmat1=matmul(umat,crcart)
        rotsph=matmul(tmat1,tmat2)
!        
!       successfull exit
!        
        return
        
      
      end subroutine rot1tosph
!EOC        
          
       
       
          
    
