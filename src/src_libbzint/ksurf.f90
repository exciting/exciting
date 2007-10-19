!BOP
!
! !ROUTINE: ksurf
!
! !INTERFACE: 
      subroutine ksurf(eo,omeg,weight)
!     
! !DESCRIPTION:
!
! This subroutine calculates the integration over an energy surface in the
! k-mesh, the area of the Fermi surface inside the tetrahedron is calculated
! which is the essential item to decide the weight of the integration over
! each vertex. 
!

! !USES:
 
! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in) :: eo(4) ! the band energies or energy difference at k
 
      real(8), intent(in) :: omeg  ! the energy surface

! !OUTPUT PARAMETERS:

      real(8), intent(out) :: weight(4) ! the weight at each corner
       
! !LOCAL VARIABLES:

      integer(4) :: i,j,k,indx
      real(8) :: zerotol, eff, denom, det, vol, total,  x1, y1
      real(8), dimension(4) :: ee,ww
      real(8), dimension(3,3) :: nodes
      real(8), dimension(3,3) :: vect 
      
! !EXTERNAL ROUTINES:
      
      real(8), external :: dos1t

! !INTRINSIC FUNCTIONS:
  
      det(i,j)=vect(2,i)*vect(3,j)-vect(2,j)*vect(3,i)

!
! !SYSTEM ROUTINES:

      intrinsic mod

! !REVISION HISTORY:
! 
! Created 17th. Jan 2005 by XZL
!
!EOP
!BOC
      zerotol=1.0d-10
      ee(1:4)=eo(1:4)
      eff=omeg
      weight(1:4)=0.0d0
      if(ee(4).le.eff) then
        indx=4
      else
        indx=0
        do i=1,4
          if(ee(i).le.eff)indx=indx+1
        enddo  
      endif
      nodes(1:3,1:3)=0.0d0
      vect(1:3,1:3)=0.0d0
      ww(1:4)=0.0d0 
      denom=(ee(2)-ee(1))**2+(ee(3)-ee(1))**2+(ee(4)-ee(1))**2
      denom=sqrt(denom)
      if(denom.lt.zerotol)goto 999
      select case(indx)
        case(0,4) 
          continue
        case(1)
!---------------------------------------------------------------------------
!   This is the case only ee(1)<efer. and we have a triangle for the surface 
!   inside tetrahedron
!--------------------------------------------------------------------------
          nodes(1,3)=(eff-ee(1))/(ee(4)-ee(1))
          nodes(2,2)=(eff-ee(1))/(ee(3)-ee(1))
          nodes(3,1)=(eff-ee(1))/(ee(2)-ee(1))
          do i=1,3
            vect(1,i)=nodes(2,i)-nodes(1,i)
            vect(2,i)=nodes(3,i)-nodes(1,i)
            vect(3,i)=(ee(i+1)-ee(1))/denom
          enddo       
          vol=0.0d0
          do i=1,3
            j=mod(i,3)+1
            k=mod(j,3)+1
            vol=vol+vect(1,i)*det(j,k)
          enddo
          vol=dabs(vol)
!--------------------------------------------------------------------------
! In above we get the coordinates of the nodes and the area of the triangle
!--------------------------------------------------------------------------

          x1=vect(2,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)- &
     &       vect(2,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          x1=x1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
          y1=0.0d0-vect(1,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)+ &
     &       vect(1,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          y1=y1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
!-------------------------------------------------------------------------
! x1, y1 corresponds to $a_1$ and $b_1$ in the introduction section about 
! integration over Fermi surface
!-------------------------------------------------------------------------        
          if(vol.eq.0.0d0)goto 999
          ww(1:4)=0.0d0
          do i=1,3
           ww(i+1)=(vect(1,i)*(1.0d0+3.0d0*x1)/6.0d0+  &
     &              vect(2,i)*(1.0d0+3.0d0*y1)/6.0d0+  &
     &              vect(3,i)*(eff-ee(1))/(2.0d0*denom))*vol
          enddo
       
          ww(1)=vol/2.0d0-ww(2)-ww(3)-ww(4)

        case(2)
!-------------------------------------------------------------------------
! This is the case when ee(1) and ee(2) < efer and we have a square for the
! surface in the tetrahedron. We further divide it into two triangles and 
! follow the same method as case(1), we get the weight on the vertices of the
! tetrahedron
!-------------------------------------------------------------------------
          nodes(1,3)=(eff-ee(1))/(ee(4)-ee(1))
          nodes(2,2)=(eff-ee(1))/(ee(3)-ee(1))
          nodes(3,1)=(ee(4)-eff)/(ee(4)-ee(2))
          nodes(3,3)=(eff-ee(2))/(ee(4)-ee(2))
      
          do i=1,3
            vect(1,i)=nodes(2,i)-nodes(1,i)
            vect(2,i)=nodes(3,i)-nodes(1,i)
            vect(3,i)=(ee(i+1)-ee(1))/denom
          enddo
       
          vol=0.0d0
          do i=1,3
            j=mod(i,3)+1
            k=mod(j,3)+1
            vol=vol+vect(1,i)*det(j,k)
          enddo
          vol=dabs(vol)
           x1=vect(2,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)- &
     &       vect(2,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          x1=x1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
          y1=0.0d0-vect(1,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)+ &
     &       vect(1,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          y1=y1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
          
          if(vol.eq.0.0d0)goto 999
          ww(1:4)=0.0d0
          do i=1,3
           ww(i+1)=(vect(1,i)*(1.0d0+3.0d0*x1)/6.0d0+  &
     &              vect(2,i)*(1.0d0+3.0d0*y1)/6.0d0+  &
     &              vect(3,i)*(eff-ee(1))/(2.0d0*denom))*vol
          enddo

          total=vol/2.0d0
     

          nodes(1:3,1:3)=0.0d0
          nodes(1,1)=(eff-ee(3))/(ee(2)-ee(3))
          nodes(1,2)=(ee(2)-eff)/(ee(2)-ee(3))
          nodes(2,2)=(eff-ee(1))/(ee(3)-ee(1))
          nodes(3,1)=(eff-ee(4))/(ee(2)-ee(4))
          nodes(3,3)=(ee(2)-eff)/(ee(2)-ee(4))

          vect(1:3,1:3)=0.0d0
          do i=1,3
            vect(1,i)=nodes(2,i)-nodes(1,i)
            vect(2,i)=nodes(3,i)-nodes(1,i)
            vect(3,i)=(ee(i+1)-ee(1))/denom
          enddo
          vol=0.0d0
          do i=1,3
            j=mod(i,3)+1
            k=mod(j,3)+1
            vol=vol+vect(1,i)*det(j,k)
          enddo
          vol=dabs(vol)
          if(vol.eq.0.0d0)goto 999
          x1=vect(2,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)- &
     &       vect(2,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          x1=x1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
          y1=0.0d0-vect(1,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)+ &
     &       vect(1,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          y1=y1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
          
        
          do i=1,3
           ww(i+1)=(vect(1,i)*(1.0d0+3.0d0*x1)/6.0d0+  &
     &              vect(2,i)*(1.0d0+3.0d0*y1)/6.0d0+  &
     &              vect(3,i)*(eff-ee(1))/(2.0d0*denom))*vol+ww(i+1)
          enddo
          total=total+vol/2.0d0
     
          ww(1)=total-ww(2)-ww(3)-ww(4)

        case(3)
!------------------------------------------------------------------------
! ee(1), ee(2) and ee(3) < efer, again, we have a triangle.
!------------------------------------------------------------------------
          nodes(1,3)=(eff-ee(1))/(ee(4)-ee(1))
          nodes(2,2)=(eff-ee(4))/(ee(3)-ee(4))
          nodes(2,3)=(ee(3)-eff)/(ee(3)-ee(4))
          nodes(3,1)=(eff-ee(4))/(ee(2)-ee(4))
          nodes(3,3)=(ee(2)-eff)/(ee(2)-ee(4))
          do i=1,3
            vect(1,i)=nodes(2,i)-nodes(1,i)
            vect(2,i)=nodes(3,i)-nodes(1,i)
            vect(3,i)=(ee(i+1)-ee(1))/denom
          enddo
          vol=0.0d0
          do i=1,3
            j=mod(i,3)+1
            k=mod(j,3)+1
            vol=vol+vect(1,i)*det(j,k)
          enddo
          vol=dabs(vol)
          if(vol.eq.0.0d0)goto 999
          x1=vect(2,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)- &
     &       vect(2,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          x1=x1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
          y1=0.0d0-vect(1,2)*(nodes(1,1)-(eff-ee(1))*vect(3,1)/denom)+ &
     &       vect(1,1)*(nodes(1,2)-(eff-ee(1))*vect(3,2)/denom)
          y1=y1/(vect(1,1)*vect(2,2)-vect(2,1)*vect(1,2))
          
          ww(1:4)=0.0d0
          do i=1,3
           ww(i+1)=(vect(1,i)*(1.0d0+3.0d0*x1)/6.0d0+  &
     &              vect(2,i)*(1.0d0+3.0d0*y1)/6.0d0+  &
     &              vect(3,i)*(eff-ee(1))/(2.0d0*denom))*vol
          enddo
          ww(1)=vol/2.0d0-ww(2)-ww(3)-ww(4)

        case default

          print*, 'case default in ksurf.f90'

      end select

      do i=1,4
        weight(i)=ww(i)/denom
      enddo

999   continue

      end subroutine ksurf

!EOC        
