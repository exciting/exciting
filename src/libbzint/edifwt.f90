
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: edifwt
!
! !INTERFACE:

      subroutine edifwt(v,omeg,figu,wt)
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

      integer(4), intent(in) :: figu ! If figu=4, it belongs to the none
!                                    equally case. If figu=6, v(1)=v(2).
!                                    If figu=8, v(1)=v(2) and v(3)=v(4).
!                                    If figu=10, v(1)=v(2)=v(3).
!                                    If figu=16, v(1)=v(2)=v(3)=v(4). 

 
! !OUTPUT PARAMETERS:            

      real(8), intent(out) :: wt  ! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: i,j

      real(8)    :: aa, bb, cc, dd, bb1, bb3,bb4
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
! Created 05.11.2004 by XZL.
!

!EOP
!BOC

      do i=1,4
       do j=1,4
        vdif(i,j)=v(i)-v(j)
       enddo
      enddo
      if(omeg.lt.1.0d-8)write(*,*)'OMEGA ZERO',omeg
      select case(figu)

      case(4)                   ! for the case none of them are equal
         
        aa=2.0d0*(omeg**2-v(1)**2)*vdif(1,2)*vdif(1,3)*vdif(1,4)*    &
     &      vdif(2,3)*vdif(2,4)*vdif(3,4) 

        dd=2.0d0*omeg*(3.0d0*v(1)**4-omeg**2*(v(3)*v(4)+v(2)*(v(3)+ &
     &     v(4)))-3.0d0*v(1)**2*(omeg**2+v(3)*v(4)+v(2)*(v(3)+v(4)))+&
     &     2.0d0*v(1)*(omeg**2*(v(3)+v(4))+v(2)*(omeg**2+3.0d0*v(3)* &
     &     v(4))))
     
        aa=aa+dd*vdif(2,3)*vdif(2,4)*vdif(3,4)*datan(v(1)/omeg)
     
        aa=aa+2.0d0*omeg*(omeg**2-3.0d0*v(2)**2)*vdif(1,3)**2*          &
     &     vdif(1,4)**2*vdif(3,4)*datan(v(2)/omeg) 
        
        aa=aa-2.0d0*omeg*(omeg**2-3.0d0*v(3)**2)*vdif(1,2)**2*          &
     &     vdif(1,4)**2*vdif(2,4)*datan(v(3)/omeg) 
     
        aa=aa+2.0d0*omeg*(omeg**2-3.0d0*v(4)**2)*vdif(1,2)**2*          &
     &     vdif(1,3)**2*vdif(2,3)*datan(v(4)/omeg) 
     
        dd=-3.0d0*omeg**2*v(2)*v(3)*v(4)+v(1)**4*(v(2)+v(3)+v(4))-      &
     &     2.0d0*v(1)**3*(3.0d0*omeg**2+v(3)*v(4)+v(2)*(v(3)+v(4)))+    &
     &     3.0d0*v(1)**2*(omeg**2*(v(3)+v(4))+v(2)*(omeg**2+v(3)*v(4)))

        bb=dd*vdif(2,3)*vdif(2,4)*vdif(3,4)*dlog(omeg**2 + v(1)**2) 
        
        bb=bb+v(2)*(3.0d0*omeg**2-v(2)**2)*vdif(1,3)**2*vdif(1,4)**2*  &
     &     vdif(3,4)*dlog(omeg**2+v(2)**2)
     
        bb=bb+v(3)*(-3.0d0*omeg**2+v(3)**2)*vdif(1,2)**2*vdif(1,4)**2*  &
     &     vdif(2,4)*dlog(omeg**2 + v(3)**2) 
     
        bb=bb-v(4)*(-3.0d0*omeg**2+v(4)**2)*vdif(1,2)**2*vdif(1,3)**2*  &
     &     vdif(2,3)*dlog(omeg**2 + v(4)**2)
     
        cc=6.0d0*vdif(1,2)**2*vdif(1,3)**2*vdif(1,4)**2*vdif(2,3)*      &
     &      vdif(2,4)*vdif(3,4)
     
      case(6)                ! for the case when v(1)=v(1)
          
        dd=v(1)**3-2.0d0*omeg**2*(v(3)+v(4))-3.0d0*v(1)**2*(v(3)+v(4))+ &
     &     v(1)*(4.0d0*omeg**2+5.0d0*v(3)*v(4))   
        
        aa=vdif(3,1)*vdif(1,4)*vdif(3,4)*dd
        
        dd=-omeg**2*(3.0d0*v(1)**2+v(3)**2+v(3)*v(4)+v(4)**2-3.0d0*v(1)*&
     &     (v(3)+v(4)))+3.0d0*(-3.0d0*v(1)**2*v(3)*v(4)+v(3)**2*v(4)**2+&
     &     v(1)**3*(v(3)+v(4)))   
        aa=aa-2.0d0*omeg*vdif(3,4)*dd*datan(v(1)/omeg)
          
        aa=aa-2.0d0*omeg*(omeg**2-3.0d0*v(3)**2)*vdif(1,4)**3*         &
     &     datan(v(3)/omeg)

        aa=aa+2.0d0*omeg*vdif(1,3)**3*(omeg**2-3.0d0*v(4)**2)*         &
     &     datan(v(4)/omeg)

        dd=-3.0d0*omeg**2*v(3)*v(4)*(v(3)+v(4))-3.0d0*v(1)**2*v(3)*    &
     &     v(4)*(v(3)+v(4))+3.0d0*v(1)*v(3)*v(4)*(3.0d0*omeg**2+v(3)*  &
     &     v(4))+v(1)**3*(-3.0d0*omeg**2+v(3)**2+v(3)*v(4)+v(4)**2)  

        bb1=-vdif(3,4)*dd
          
        bb3=v(3)*vdif(1,4)**3*(v(3)**2-3.0d0*omeg**2)
     
        bb4=-v(4)*vdif(1,3)**3*(v(4)**2-3.0d0*omeg**2)
 
        bb=bb1*dlog((v(1)**2+omeg**2))+               &
     &     bb3*dlog((v(3)**2+omeg**2))+               &
     &     bb4*dlog(v(4)**2+omeg**2)
      
        cc=6.0d0*vdif(1,3)**3*vdif(1,4)**3*vdif(3,4)

      case(8)          !for the case when v(1)=v(1) and v(3)=v(4)
         vp=(v(3)+v(4))*0.5d0
         ev(1)=v(1)
         ev(2)=v(2)
         ev(3)=vp
         ev(4)=vp
         do i=1,4
           do j=1,4
             vdif(i,j)=ev(i)-ev(j)
           enddo
         enddo
        dd=6.0d0*omeg**2+ev(1)**2-5.0d0*ev(1)*ev(3)-2.0d0*ev(3)**2
        aa=vdif(3,1)*dd
  
        dd=6.0d0*omeg*(omeg**2-ev(3)*(ev(3)+2.0d0*ev(1)))
        aa=aa+dd*(datan(ev(1)/omeg)-datan(ev(3)/omeg))

            bb=-3.0d0*(ev(1)*ev(3)**2-(ev(1)+2.0d0*ev(3))*omeg**2)
            bb=bb*dlog((ev(1)**2+omeg**2)/(ev(3)**2+omeg**2))

         
        cc=6.0d0*vdif(1,3)**4

      case(10)             ! for the case when v(1)=v(1)=v(3)

        aa=vdif(1,4)*(6.0d0*omeg**2-2.0d0*v(1)**2+7.0d0*v(1)*v(4)- &
     &     11.0d0*v(4)**2)     
          
        dd=6.0d0*omeg*(omeg**2-3.0d0*v(4)**2)
        aa=aa-dd*(datan(v(1)/omeg)-datan(v(4)/omeg))
          
        dd=3.0d0*v(4)*(v(4)**2-3.0d0*omeg**2)
        bb=dd*dlog((v(1)**2+omeg**2)/(v(4)**2+omeg**2))
         
        cc=18.0d0*vdif(1,4)**4

      case(16)
        aa=-v(1)
        bb=0.0d0
        cc=12.0d0*(omeg**2+v(1)**2)

      end select
!      wt=0.0d0
!      if((dabs(cc)*1.0d+2).gt.dabs(aa+bb)) then
       wt=(aa+bb)/cc
!      endif
      if(abs(wt).gt.1.0d+1)then
        write(*,1)wt,figu,omeg,v,vdif,aa,bb,cc
!        stop
      endif
    1 format('warning: weightt =',g18.10,' case:',i4,/,'omeg =',g18.10,/,'v =',4g18.10,    &
     &    /,'vdif = ',/,4(4g18.10,/),' aa =',g18.10,' bb =',g18.10,      &
     &    ' cc =',g18.10)     

      return
      
      end subroutine edifwt
!EOC      
