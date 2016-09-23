!BOP
!
! !ROUTINE: ksurf
!
! !INTERFACE: 
      subroutine ksurf(ei,esurf,weight)
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

      real(8), intent(in) :: ei(4) ! the band energies or energy difference at k
 
      real(8), intent(in) :: esurf  ! the energy or energy differnce surface

! !OUTPUT PARAMETERS:

      real(8), intent(out) :: weight(4) ! the weight at each corner
       
! !LOCAL VARIABLES:
     
     integer(4) :: i, j
     integer(4) :: indx ! Number of nodes below esurf
     integer(4) :: sort_ind(4)
            
     real(8) :: delta21, delta31, delta41, delta32, delta42, delta43  ! energy differences       
     real(8) :: deles1, deles2, deles3, deles4 ! esurf-eo(i) 
     real(8) :: den0, den1, den02, den12 ! denominators
     real(8) :: num0, num1, num02, num12 ! numerators
     real(8) :: wt         ! total weight
     real(8) :: eo(4), wtmp(4)
     
! !EXTERNAL ROUTINES:
      

! !INTRINSIC FUNCTIONS:
  

!
! !SYSTEM ROUTINES:

! !REVISION HISTORY:
! 
! Created 17th. Jan 2005 by XZL
! Last version: May 2007 by RGA
!EOP
!BOC
      wtmp(1:4)=0.0d0
      sort_ind(1)=minloc(ei,dim=1)
      sort_ind(4)=maxloc(ei,dim=1)
      if(sort_ind(1).ne.sort_ind(4))then
        j=2
        do i=1,4
          if((i.ne.sort_ind(1)).and.(i.ne.sort_ind(4)))then
            sort_ind(j)=i
            j=j+1
          endif
        enddo
        if(ei(sort_ind(2)).gt.ei(sort_ind(3)))then
          j=sort_ind(2)
          sort_ind(2)=sort_ind(3)
          sort_ind(3)=j
        endif
        do i=1,4  
          eo(i)=ei(sort_ind(i))
        enddo
      else
        eo=ei
        do i=1,4
          sort_ind(i)=i
        enddo      
      endif  
        
      if(esurf.lt.eo(1))then
        indx=0
      else
        indx=1
        do i=2,4
          if(eo(i).le.esurf)indx=indx+1
        enddo  
      endif
      select case(indx)
        case(0,4) 
          continue
        case(1) ! Only eo(1)< esurf, one triangle

          deles1=esurf-eo(1)
          delta21=eo(2)-eo(1)
          delta31=eo(3)-eo(1)
          delta41=eo(4)-eo(1)
          
          den0=6.0d0*delta21*delta31*delta41
          num0=deles1*deles1
          
! total wtmp
          wt= 3.0d0*num0/den0
          num1=num0*deles1
! wtmp(2)          
          den1=delta21*den0
          wtmp(2)=num1/den1
! wtmp(3)
          den1=delta31*den0
          wtmp(3)=num1/den1
! wtmp(4)
          den1=delta41*den0
          wtmp(4)=num1/den1
! wtmp(1)          
          wtmp(1)=wt-wtmp(2)-wtmp(3)-wtmp(4)

        case(2) ! eo(1)<eo(2)<esurf, two triangles
        
          deles1=esurf-eo(1)
          deles2=esurf-eo(2)
          deles3=eo(3)-esurf
          deles4=eo(4)-esurf
          delta21=eo(2)-eo(1)
          delta31=eo(3)-eo(1)
          delta41=eo(4)-eo(1)
          delta32=eo(3)-eo(2)
          delta42=eo(4)-eo(2)
          delta43=eo(4)-eo(3)
          
! total wtmp          
          den0=6.0d0*delta31*delta32*delta41
          den1=6.0d0*delta32*delta41*delta42
          
          num0=deles1*deles3
          num1=deles2*deles4
          
          wt=3.0d0*(num0/den0+num1/den1)
          
! wtmp(2)
         den02=delta32*den0
         den12=delta41*delta42*den1
         num02=deles3*num0
         num12=(deles3*delta42+deles2*delta31)*num1

         wtmp(2)=num02/den02+num12/den12

! wtmp(3)
         den02=delta31*den02
         den12=delta32*den1
         num02=(deles1*delta32+deles2*delta31)*num0
         num12=deles2*num1

         wtmp(3)=num02/den02+num12/den12

! wtmp(4)
         den02=delta41*den0
         den12=delta41*delta42*den1
         num02=deles1*num0
         num12=(deles1*delta42+deles2*delta41)*num1

         wtmp(4)=num02/den02+num12/den12

! wtmp(1)
          wtmp(1)=wt-wtmp(2)-wtmp(3)-wtmp(4)

        case(3) ! ee(1), ee(2) and ee(3) < efer, one triangle.

          deles4=eo(4)-esurf
          delta41=eo(4)-eo(1)
          delta42=eo(4)-eo(2)
          delta43=eo(4)-eo(3)
          
          den0=6.0d0*delta41*delta42*delta43
          if(den0.lt.1.0e-12) then
            write(6,*)'ERROR in ksurf: flat band'
            write(6,*) eo
            stop 'flat band'
          endif    
          num0=deles4*deles4
          
! total wtmp
          wt= 3.0d0*num0/den0
          num1=num0*deles4
! wtmp(1)          
          den1=delta41*den0
          wtmp(1)=num1/den1
! wtmp(2)
          den1=delta42*den0
          wtmp(2)=num1/den1
! wtmp(3)
          den1=delta43*den0
          wtmp(3)=num1/den1
! wtmp(4)          
          wtmp(4)=wt-wtmp(1)-wtmp(2)-wtmp(3)

      end select
      do i=1,4
        weight(sort_ind(i))=wtmp(i)
      enddo  

      end subroutine ksurf

!EOC        
