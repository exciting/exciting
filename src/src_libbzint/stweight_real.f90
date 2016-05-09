!BOP
!
! !ROUTINE: redifwt
!
! !INTERFACE:

      subroutine stweight_real(deltae_vert,freq,equiv_flag,weight_vert)

! !DESCRIPTION:
!
! This subroutine calculates the weight on the second vertex of the tetrahedron 
! which is full occupied for k and full unoccupied for k-q. It is also for the 
! case of $sigfreq=2$ when we consider the real frequency contribution.
! 
                                            
! !USES:
      use tetra_internal, only: weighttol,weightwarn,fout

! !INPUT PARAMETERS:
      implicit none
      
      real(8), intent(in) :: deltae_vert(4)  ! energy differences at four vertices 
      real(8), intent(in) :: freq            ! the frequency omega to be calculated

      integer(4), intent(in) :: equiv_flag   ! equiv_flag = 4, none is equal 
                                             ! equiv_flag=6, v(1)=v(2).
                                             ! equiv_flag=8, v(1)=v(2) and v(3)=v(4).
                                             ! equiv_flag=10, v(1)=v(2)=v(3).
                                             ! e1uiv_flag=16, v(1)=v(2)=v(3)=v(4). 

 
! !OUTPUT PARAMETERS:            

      real(8), intent(out) :: weight_vert(4)      ! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: i,j, k
      integer(4) :: ivert
      integer(4) :: isign_om

      real(8)    :: aa,  cc, dd, vp
      real(8)    :: omeg
      
      real(8), dimension(4,4) :: vdif

      real(8), dimension(4) :: ev
      real(8), dimension(2,4) :: weight_tmp
      logical:: lwarn = .False.
      character(20):: sname="stweight_real"

! !DEFINED PARAMETERS:

      real(8), parameter :: haier=1.0d-12
! 
! !SYSTEM ROUTINES:

      intrinsic dlog
      intrinsic dabs
 
! !REVISION HISTORY:
!
! Created 05.11.2004 by XZL.
!
! Intrinsic functions

!EOP
!BOC
      weight_tmp=0.0d0
      
      do isign_om=1,2
        omeg=dble(3-2*isign_om)*freq

        select case(equiv_flag)

        case(4)                        ! for the case of none of them are equal

          do ivert=1,4
            do j=1,4
              k=mod(j+ivert-2,4)+1       
              ev(k)=deltae_vert(j)
            enddo
            call set_vdif
            
            if(dabs(omeg-ev(1)).lt.haier) then
 
              aa=-dlog(dabs(vdif(2,1)))*vdif(2,1)*vdif(4,3)

              aa=aa+dlog(dabs(vdif(3,1)))*vdif(3,1)*vdif(4,2)

              aa=aa-dlog(dabs(vdif(4,1)))*vdif(4,1)*vdif(3,2)

              cc=6.0d0*vdif(3,2)*vdif(4,3)*vdif(4,2)

            else if(dabs(omeg-ev(2)).lt.haier) then
         
              aa=dlog(dabs(vdif(3,2)))*vdif(3,2)**2*vdif(4,1)**2

              aa=aa-dlog(dabs(vdif(4,2)))*vdif(3,1)**2*vdif(4,2)**2

              dd=dlog(dabs(vdif(2,1)))*(vdif(3,1)*vdif(4,2)+vdif(3,2)*    &
     &           vdif(4,1))
              dd=dd+vdif(3,1)*vdif(4,1)
 
              aa=aa+vdif(2,1)*vdif(4,3)*dd

              cc=6.0d0*vdif(3,1)**2*vdif(4,1)**2*vdif(4,3)

            else if(dabs(omeg-ev(3)).lt.haier) then
    
              aa=dlog(dabs(vdif(3,2)))*vdif(3,2)**2*vdif(4,1)**2

              aa=aa-dlog(dabs(vdif(4,3)))*vdif(2,1)**2*vdif(4,3)**2
 
              dd=dlog(dabs(vdif(3,1)))*(vdif(4,1)*vdif(3,2)-vdif(4,3)*  &
     &           vdif(2,1))
              dd=dd-vdif(2,1)*vdif(4,1)

              aa=aa-vdif(3,1)*vdif(4,2)*dd

              cc=6.0d0*vdif(2,1)**2*vdif(4,1)**2*vdif(4,2)
  
            else if(dabs(omeg-ev(4)).lt.haier) then

              aa=0.0d0-dlog(dabs(vdif(4,3)))*vdif(3,1)**2*vdif(4,2)**2

              aa=aa+dlog(dabs(vdif(4,2)))*vdif(4,2)**2*vdif(3,1)**2

              dd=dlog(dabs(vdif(4,1)))-dlog(dabs(vdif(4,3)))
              dd=dd*(vdif(3,1)*vdif(4,2)+vdif(2,1)*vdif(4,3))-          &
     &           vdif(2,1)*vdif(3,1)

              aa=aa-vdif(3,2)*vdif(4,1)*dd

              cc=6.0d0*vdif(2,1)**2*vdif(3,1)**2*vdif(3,2)

            else

              dd=(omeg-ev(4))**3*vdif(2,1)**2*vdif(3,2)*vdif(3,1)**2
    
              aa=dlog(dabs(omeg-ev(4)))*dd

              dd=(omeg-ev(3))**3*vdif(2,1)**2*vdif(4,2)*vdif(4,1)**2

              aa=aa-dlog(dabs(omeg-ev(3)))*dd
  
              dd=-vdif(2,1)*vdif(3,1)*(omeg-ev(4))-(omeg-ev(2))*vdif(3,1)*   &
     &           vdif(4,1)
              dd=dd-(omeg-ev(3))*vdif(2,1)*vdif(4,1)
              dd=dd*vdif(3,2)*vdif(4,2)*vdif(4,3)*(omeg-ev(1))**2

              aa=aa+dlog(dabs(omeg-ev(1)))*dd

              aa=aa+(omeg-ev(1))**2*vdif(2,1)*vdif(3,2)*vdif(3,1)*         &
     &           vdif(4,2)*vdif(4,1)*vdif(4,3)

              aa=aa+dlog(dabs(omeg-ev(2)))*(omeg-ev(2))**3*vdif(3,1)**2*    &
     &           vdif(4,1)**2*vdif(4,3)

              cc=6.0d0*vdif(2,1)**2*vdif(3,2)*vdif(3,1)**2*vdif(4,2)*     &
     &           vdif(4,1)**2*vdif(4,3)
          
            endif
            call sub_set_weight
          enddo ! ivert  

        case(6)                     ! for the case that a=b
  
          do ivert=1,2
            ev(1)=deltae_vert(ivert)
            ev(2)=ev(1)
            ev(3:4)=deltae_vert(3:4)
            call set_vdif
    
            if(dabs(omeg-ev(1)).lt.haier) then
  
              aa=dlog(dabs(vdif(3,1)))-dlog(dabs(vdif(4,1)))
   
              cc=6.0d0*vdif(4,3)
  
            else if(dabs(omeg-ev(3)).lt.haier) then
   
              aa=2.0d0*(dlog(dabs(vdif(3,1)))-dlog(dabs(vdif(4,1))))*     &
     &           vdif(4,1)**2

              aa=aa+(2.0d0*vdif(4,3)+vdif(4,1))*vdif(4,1)

              cc=12.0d0*vdif(4,1)**3

            else if(dabs(omeg-ev(4)).lt.haier) then
 
              aa=vdif(3,1)*(vdif(3,1)-2.0d0*vdif(4,3))
              aa=aa-2.0d0*(dlog(dabs(vdif(4,3)))-dlog(dabs(vdif(4,1))))*  &
     &           vdif(4,3)**2

              cc=12.0d0*vdif(3,1)**3

            else

              aa=2.0d0*dlog(dabs(omeg-ev(4)))*(omeg-ev(4))**3*            &
     &           vdif(3,1)**3
        
              aa=aa-2.0d0*dlog(dabs(omeg-ev(3)))*(omeg-ev(3))**3*         &
     &           vdif(4,1)**3

              dd=(omeg-ev(4))**2*vdif(3,1)**2+(omeg-ev(3))*vdif(3,1)*     &
     &           (omeg-ev(4))*vdif(4,1)
              dd=dd+(omeg-ev(3))**2*vdif(4,1)**2

              aa=aa+2.0d0*dlog(dabs(omeg-ev(1)))*(omeg-ev(1))*dd*vdif(4,3)

              dd=vdif(3,1)*vdif(4,1)-2.0d0*vdif(3,1)*(omeg-ev(4))-2.0d0*  &
     &           (omeg-ev(3))*vdif(4,1)

              aa=aa+(omeg-ev(1))*dd*vdif(3,1)*vdif(4,1)*vdif(4,3)

              cc=12.0d0*vdif(3,1)**3*vdif(4,1)**3*vdif(4,3)
            
            endif
            call sub_set_weight
          enddo ! ivert  

          do ivert=3,4
            ev(1)=(deltae_vert(1)+deltae_vert(2))*0.5d0
            ev(2)=ev(1)
            do j=3,4
              k=mod(j+ivert,2)+3
              ev(k)=deltae_vert(j)
            enddo
            call set_vdif
  
            if(dabs(omeg-ev(1)).lt.haier) then

              aa=(dlog(dabs(vdif(3,1)))-dlog(dabs(vdif(4,1))))*vdif(4,1)
              aa=aa+vdif(4,3)
 
              cc=6.0d0*vdif(4,3)**2

            else if(dabs(omeg-ev(3)).lt.haier) then
 
              aa=(dlog(dabs(vdif(3,1)))-dlog(dabs(vdif(4,3))))*vdif(4,3)
              aa=aa+vdif(4,1)

              cc=6.0d0*vdif(4,1)**2

            else if(dabs(omeg-ev(4)).lt.haier) then
 
              aa=vdif(3,1)*(vdif(4,3)+vdif(4,1))
              aa=aa-2.0d0*(dlog(dabs(vdif(4,1)))-dlog(dabs(vdif(4,3))))*  &
     &           vdif(4,1)*vdif(4,3)

              cc=6.0d0*vdif(3,1)**3

            else

              aa=dlog(dabs(omeg-ev(4)))*(omeg-ev(4))**3*vdif(3,1)**3
                    
              dd=2.0d0*vdif(4,3)*(omeg-ev(1))-(omeg-ev(4))*vdif(3,1)
              aa=aa+dlog(dabs(omeg-ev(3)))*(omeg-ev(3))**2*dd*vdif(4,1)**2

              dd=(omeg-ev(3))**2*vdif(4,1)+(omeg-ev(1))**2*vdif(4,3)
              aa=aa+dd*vdif(3,1)*vdif(4,1)*vdif(4,3)

              dd=2.0d0*(omeg-ev(3))*vdif(4,1)+vdif(3,1)*(omeg-ev(4))
              aa=aa-dlog(dabs(omeg-ev(1)))*(omeg-ev(1))**2*dd*vdif(4,3)**2

              cc=6.0d0*vdif(3,1)**3*vdif(4,1)**2*vdif(4,3)**2
            
            endif
            call sub_set_weight
          enddo ! ivert  
            
        case(8)                      ! for the case that a=b and c=d

          do ivert=1,2
            ev(1)=deltae_vert(ivert)
            ev(2)=ev(1)
            vp=(deltae_vert(3)+deltae_vert(4))*0.5d0
            ev(3)=vp
            ev(4)=vp
            call set_vdif
            
            if(dabs(omeg-ev(1)).lt.haier) then
          
              aa=-1.0d0
              cc=6.0d0*vdif(3,1)

            else if(dabs(omeg-ev(3)).lt.haier) then

              aa=1.0d0
              cc=12.0d0*vdif(3,1)

            else

              dd=dlog(dabs(omeg-ev(1)))-dlog(dabs(omeg-ev(3)))

              aa=6.0d0*(omeg-ev(1))*dd*(omeg-ev(3))**2
           
              dd=9.0d0*(omeg-ev(1))*vdif(3,1)-2.0d0*vdif(3,1)**2-6.0d0*   &
     &           (omeg-ev(1))**2
            
              aa=aa+vdif(3,1)*dd
        
              cc=12.0d0*vdif(3,1)**4    

            endif
            call sub_set_weight
          enddo ! ivert  

          do ivert=3,4
            ev(3)=deltae_vert(ivert)
            ev(4)=ev(3)
            vp=(deltae_vert(1)+deltae_vert(2))*0.5d0
            ev(1)=vp
            ev(2)=vp
            call set_vdif

            if(dabs(omeg-ev(1)).lt.haier) then
          
              aa=-1.0d0
            
              cc=12.0d0*vdif(3,1)

            else if(dabs(omeg-ev(3)).lt.haier) then

              aa=1.0d0
            
              cc=6.0d0*vdif(3,1)
 
            else

              dd=dlog(dabs(omeg-ev(3)))-dlog(dabs(omeg-ev(1)))

              aa=6.0d0*(omeg-ev(1))**2*dd*(omeg-ev(3))
           
              dd=3.0d0*(omeg-ev(1))*vdif(3,1)+vdif(3,1)**2-6.0d0*          &
     &           (omeg-ev(1))**2
              aa=aa-vdif(3,1)*dd
        
              cc=12.0d0*vdif(3,1)**4    

            endif
            call sub_set_weight
          enddo ! ivert  
            
        case(10)                   ! for the case that a=b=c

          do ivert=1,3
            ev(1)=deltae_vert(ivert)
            ev(2)=ev(1)
            ev(3)=ev(1)
            ev(4)=deltae_vert(4)
            call set_vdif

            if(dabs(omeg-ev(4)).lt.haier) then
              aa=1.0d0
              cc=18.0d0*vdif(4,1)
            else
              aa=6.0d0*(dlog(dabs(omeg-ev(4)))-dlog(dabs(omeg-ev(1))))* &
     &           (omeg-ev(4))**3
              dd=0.0d0-15.0d0*(omeg-ev(1))*vdif(4,1)+11.0d0*vdif(4,1)**2+ &
     &           6.0d0*(omeg-ev(1))**2
              aa=aa+vdif(4,1)*dd
              cc=36.0d0*vdif(4,1)**4
            endif
            call sub_set_weight
          enddo ! ivert  

          vp=sum(deltae_vert(1:3))/3.0d0
          ev(1:3)=vp
          ev(4)=deltae_vert(4) 
          call set_vdif
          ivert = 4
          if(dabs(omeg-ev(4)).lt.haier) then
            aa=1.0d0
            cc=12.0d0*vdif(4,1)
          else
            dd=5.0d0*ev(1)*ev(4)+2.0d0*ev(4)**2-3.0d0*ev(1)*omeg-ev(1)**2
            dd=dd-9.0d0*ev(4)*omeg+6.0d0*omeg**2
            aa=vdif(1,4)*dd
  
            dd=dlog(dabs(omeg-ev(1)))-dlog(dabs(omeg-ev(4)))
            aa=aa-6.0d0*(ev(1)-omeg)*(ev(4)-omeg)**2*dd
 
            cc=12.0d0*vdif(4,1)**4
          endif 
          call sub_set_weight

        case(16)                   ! for the case that a=b=c=d
          ev(1:4)=deltae_vert(1:4)
          aa=1.0d0
          do ivert=1,4    
            cc=24.0d0*(omeg-ev(ivert))
            call sub_set_weight
          enddo ! ivert
        end select
        
      enddo ! isign_om
      do ivert=1,4
        weight_vert(ivert)=weight_tmp(1,ivert)+weight_tmp(2,ivert)
      enddo  

      contains
        subroutine sub_set_weight
          if((dabs(cc)*weighttol).gt.dabs(aa)) then
            weight_tmp(isign_om,ivert)=aa/cc
          else
            if(abs(cc)<1.e-20) cc = 1.e-20
            weight_tmp(isign_om,ivert)=aa/cc
            write(fout,1)ivert,weight_tmp(isign_om,ivert),equiv_flag,&
     &                 omeg,deltae_vert,vdif,aa,cc
          endif

    1     format('warning in stweight_real: large weight !!!',/, &
     &    'vertix:',i4,' weightt =',g18.10,&
     &    ' case:',i4,/,'omeg =',g18.10,/,'deltae_vert =',4g18.10,/,    &
     &    'vdif = ',/,4(4g18.10,/),' aa =',g18.10,                      &
     &    ' cc =',g18.10)     
        end subroutine 
      
        subroutine set_vdif
           do i=1,4
            do j=1,4
              vdif(i,j)=ev(i)-ev(j)
            enddo
          enddo ! i  
        end subroutine set_vdif
      
      end subroutine stweight_real
!EOC
