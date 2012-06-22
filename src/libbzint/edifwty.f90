
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: edifwty
!
! !INTERFACE:

      subroutine edifwty(v,omeg,figu,wt)

! !DESCRIPTION:
!
! This subroutine calculates the weight on vertex 3 of the 
! small tetrahedron. For the case of $sigfreq=3$ when we use 
! the imaginary frequency.
!  
                                            
! !INPUT PARAMETERS:
      implicit none
      
      real(8), intent(in) :: v(4)   ! difference of the energy
!                           in k-mesh tetrahedron vertices and 
!                           k-q mesh tetrahedron vertices.

      real(8), intent(in) :: omeg  ! the frequency omega to be
!                                    calculated

      integer(4), intent(in) :: figu  ! If figu=4, it belongs 
!                          to the none equally case. If figu=6,
!                          v(1)=v(2). If figu=8, v(1)=v(2) and 
!                          v(3)=v(4). If figu=10, v(1)=v(2)=v(3).
!                          If figu=16, v(1)=v(2)=v(3)=v(4). 

 
! !OUTPUT PARAMETERS:            

      real(8), intent(out) :: wt   ! the weight on vertex 4.

! !LOCAL VARIABLES:

      integer(4) :: i,j

      real(8)    :: aa, bb, cc, dd
      real(8)    :: vp 

      real(8), dimension(4) :: ev
      real(8), dimension(4,4) :: vdif

! !DEFINED PARAMETERS:
   
      real(8), parameter :: haier=1.0d-20

! 
! !SYSTEM ROUTINES:

   intrinsic datan
   intrinsic dlog

 
! !REVISION HISTORY:
!
! Created 04.11.2004 by XZL.
!
!EOP
!BOC


       select case(figu)

       case(4)    ! for the case when none of them are equal
         do i=1,4
           do j=1,4
             vdif(i,j)=v(i)-v(j)
           enddo
         enddo
         
         aa=2.0d0*(v(3)**2-omeg**2)*vdif(1,2)*vdif(1,3)*vdif(2,3)*     &
     &      vdif(1,4)*vdif(2,4)*vdif(3,4)
          
         aa=aa-2.0d0*omeg*(omeg**2-3.0d0*v(1)**2)*vdif(2,3)**2*        &
     &      vdif(2,4)*vdif(3,4)**2*datan(v(1)/omeg)

         aa=aa+2.0d0*omeg*(omeg**2-3.0d0*v(2)**2)*vdif(1,3)**2*        &
     &      vdif(1,4)*vdif(3,4)**2*datan(v(2)/omeg)

         dd=omeg**2*(-2.0d0*v(2)*v(3)+3.0d0*v(3)**2+v(2)*v(4)-     &
     &      2.0d0*v(3)*v(4)+v(1)*(v(2)-2.0d0*v(3)+v(4)))    
         dd=dd+3.0d0*v(3)*(-v(3)**3+v(2)*v(3)*v(4)+v(1)*(v(2)*   &
     &      v(3)-2.0d0*v(2)*v(4)+v(3)*v(4)))    
        
         aa=aa+2.0d0*omeg*vdif(1,2)*vdif(1,4)*vdif(2,4)*dd*             &
     &      datan(v(3)/omeg)
          
         aa=aa-2.0d0*omeg*vdif(1,2)*vdif(1,3)**2*vdif(2,3)**2*          &  
     &      (omeg**2-3.0d0*v(4)**2)*datan(v(4)/omeg)

         bb=v(1)*(v(1)**2-3.0d0*omeg**2)*vdif(2,3)**2*vdif(2,4)*      &
     &      vdif(3,4)**2*dlog(v(1)**2+omeg**2)
     
         bb=bb-v(2)*(v(2)**2-3*omeg**2)*vdif(1,3)**2*vdif(1,4)*       &
     &      vdif(3,4)**2*dlog(v(2)**2+omeg**2)
     
         dd=-3.0d0*omeg**2*(v(3)**2*(v(2)-2.0d0*v(3)+v(4))+v(1)* &
     &      (v(3)**2-v(2)*v(4)))
         dd=dd-v(3)**2*(v(1)*(-2.0d0*v(2)*v(3)+v(3)**2+3.0d0*  &
     &      v(2)*v(4)-2.0d0*v(3)*v(4))+v(3)*(v(2)*v(3)-2.0d0*    &
     &      v(2)*v(4)+v(3)*v(4)))    

         bb=bb+vdif(1,2)*vdif(1,4)*vdif(2,4)*dd*dlog(v(3)**2+omeg**2)

         bb=bb+v(4)*vdif(1,2)*vdif(1,3)**2*vdif(3,2)**2*         &
     &      (v(4)**2-3.0d0*omeg**2)*dlog(v(4)**2+omeg**2)
         
         cc=6.0d0*vdif(1,2)*vdif(1,3)**2*vdif(2,3)**2*vdif(1,4)* &
     &      vdif(4,2)*vdif(3,4)**2
    
       case(6)               ! for the case when v(1)=v(2)
         vp=(v(1)+v(2))*0.5d0
         ev(1)=vp
         ev(2)=vp
         ev(3)=v(3)
         ev(4)=v(4)
         do i=1,4
           do j=1,4
             vdif(i,j)=ev(i)-ev(j)
           enddo
         enddo
          
         dd=-ev(1)*ev(3)**2+omeg**2*(ev(1)+ev(3)-2.0d0*ev(4))+ev(3)**2* & 
     &      ev(4)+ev(1)**2*vdif(4,3)
         aa=2.0d0*dd*vdif(1,3)*vdif(1,4)*vdif(3,4)
          
         dd=3.0d0*ev(1)**3+3.0d0*ev(1)**2*ev(3)+omeg**2*(ev(3)+2.0d0*   &
     &      ev(4))-3.0d0*ev(1)*(omeg**2+2.0d0*ev(3)*ev(4))
         aa=aa+2.0d0*vdif(3,4)**2*omeg*dd*datan(ev(1)/omeg)

         dd=3.0d0*ev(3)*(ev(3)**2+ev(1)*(ev(3)-2.0d0*ev(4)))+omeg**2*   &
     &      (ev(1)-3.0d0*ev(3)+2.0d0*ev(4))    
         aa=aa-2.0d0*vdif(1,4)**2*omeg*dd*datan(ev(3)/omeg)

         aa=aa+2.0d0*vdif(1,3)**3*omeg*(omeg**2-3.0d0*ev(4)**2)*        &
     &      datan(ev(4)/omeg)

         dd=ev(1)**2*(-3.0d0*ev(3)*ev(4)+ev(1)*(2.0d0*ev(3)+ev(4)))
         dd=dd+omeg**2*(-6.0d0*ev(1)**2+3.0d0*ev(1)*ev(4)+3.0d0*ev(3)*  &
     &      ev(4))
         bb=vdif(3,4)**2*dd*dlog(ev(1)**2+omeg**2)

         dd=ev(3)**2*(2.0d0*ev(1)*ev(3)-3.0d0*ev(1)*ev(4)+ev(3)*ev(4))
         dd=dd+3.0d0*omeg**2*(-2.0d0*ev(3)**2+(ev(1)+ev(3))*ev(4))
        
         bb=bb-vdif(1,4)**2*dd*dlog(ev(3)**2+omeg**2)

         bb=bb-ev(4)*vdif(1,3)**3*(ev(4)**2-3.0d0*omeg**2)*             &
     &      dlog(ev(4)**2+omeg**2)
         
         cc=6.0d0*vdif(1,3)**3*vdif(1,4)**2*vdif(3,4)**2

      case(8)      ! for the case when v(1)=v(2) and v(3)=v(4)
         vp=(v(1)+v(2))*0.5d0
         ev(1)=vp
         ev(2)=vp
         ev(3)=v(3)
         ev(4)=v(4)
         do i=1,4
           do j=1,4
             vdif(i,j)=ev(i)-ev(j)
           enddo
         enddo
          
          dd=6.0d0*omeg**2-2.0d0*ev(1)**2-5.0d0*ev(1)*ev(3)+ev(3)**2
          aa=vdif(1,3)*dd
          
          dd=6.0d0*omeg*(ev(1)**2+2.0d0*ev(1)*ev(3)-omeg**2)
          aa=aa+dd*(datan(ev(1)/omeg)-datan(ev(3)/omeg))

          dd=3.0d0*(ev(1)**2*ev(3)-(2.0d0*ev(1)+ev(3))*omeg**2)
          bb=dd*dlog((ev(1)**2+omeg**2)/(ev(3)**2+omeg**2))

         
          cc=6.0d0*vdif(1,3)**4

      case(10)         ! for the case when v(1)=v(2)=v(3)
         do i=1,4
           do j=1,4
             vdif(i,j)=v(i)-v(j)
           enddo
         enddo
         
          aa=vdif(3,4)*(6.0d0*omeg**2-2.0d0*v(3)**2+7.0d0*v(3)*v(4)- &
     &       11.0d0*v(4)**2)     
          
          dd=6.0d0*omeg*(omeg**2-3.0d0*v(4)**2)
          aa=aa-dd*(datan(v(3)/omeg)-datan(v(4)/omeg))
          
          dd=3.0d0*v(4)*(v(4)**2-3.0d0*omeg**2)
          bb=dd*dlog((v(3)**2+omeg**2)/(v(4)**2+omeg**2))

         
          cc=18.0d0*vdif(3,4)**4

       case(16)
        aa=-v(3)
        bb=0.0d0
        cc=12.0d0*(omeg**2+v(3)**2)
      end select
!      wt=0.0d0
!      if((dabs(cc)*1.0d+2).gt.dabs(aa+bb)) then
        wt=(aa+bb)/cc
!      endif
      if(abs(wt).gt.1.0d+1)then
        write(*,1)wt,figu,omeg,v,vdif,aa,bb,cc
      endif
    1 format('warning: weighty =',g18.10,' case:',i4,/,'omeg =',g18.10,/,'v =',4g18.10,    &
     &    /,'vdif = ',/,4(4g18.10,/),' aa =',g18.10,' bb =',g18.10,      &
     &    ' cc =',g18.10)     

      return
      
      end subroutine edifwty
!EOC      
