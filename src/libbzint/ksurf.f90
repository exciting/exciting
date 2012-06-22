
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: ksurf
!
! !INTERFACE: 
      subroutine ksurf(eo,esurf,weight)
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
 
      real(8), intent(in) :: esurf  ! the energy or energy differnce surface

! !OUTPUT PARAMETERS:

      real(8), intent(out) :: weight(4) ! the weight at each corner
       
! !LOCAL VARIABLES:
     
     integer(4) :: i 
     integer(4) :: indx ! Number of nodes below esurf
            
     real(8) :: delta21, delta31, delta41, delta32, delta42, delta43  ! energy differences       
     real(8) :: deles1, deles2, deles3, deles4 ! esurf-eo(i) 
     real(8) :: den0, den1, den02, den12 ! denominators
     real(8) :: num0, num1, num02, num12 ! numerators
     real(8) :: wt         ! total weight
     
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
      weight(1:4)=0.0d0
      
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
          
! total weight
          wt= 3.0d0*num0/den0
          num1=num0*deles1
! weight(2)          
          den1=delta21*den0
          weight(2)=num1/den1
! weight(3)
          den1=delta31*den0
          weight(3)=num1/den1
! weight(4)
          den1=delta41*den0
          weight(4)=num1/den1
! weight(1)          
          weight(1)=wt-weight(2)-weight(3)-weight(4)

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
          
! total weight          
          den0=6.0d0*delta31*delta32*delta41
          den1=6.0d0*delta32*delta41*delta42
          
          num0=deles1*deles3
          num1=deles2*deles4
          
          wt=3.0d0*(num0/den0+num1/den1)
          
! weight(2)
         den02=delta32*den0
         den12=delta41*delta42*den1
         num02=deles3*num0
         num12=(deles3*delta42+deles2*delta31)*num1

         weight(2)=num02/den02+num12/den12

! weight(3)
         den02=delta31*den02
         den12=delta32*den1
         num02=(deles1*delta32+deles2*delta31)*num0
         num12=deles2*num1

         weight(3)=num02/den02+num12/den12

! weight(4)
         den02=delta41*den0
         den12=delta41*delta42*den1
         num02=deles1*num0
         num12=(deles1*delta42+deles2*delta41)*num1

         weight(4)=num02/den02+num12/den12

! weight(1)
          weight(1)=wt-weight(2)-weight(3)-weight(4)

        case(3) ! ee(1), ee(2) and ee(3) < efer, one triangle.

          deles4=eo(4)-esurf
          delta41=eo(4)-eo(1)
          delta42=eo(4)-eo(2)
          delta43=eo(4)-eo(3)
          
          den0=6.0d0*delta41*delta42*delta43
          num0=deles4*deles4
          
! total weight
          wt= 3.0d0*num0/den0
          num1=num0*deles4
! weight(1)          
          den1=delta41*den0
          weight(1)=num1/den1
! weight(2)
          den1=delta42*den0
          weight(2)=num1/den1
! weight(3)
          den1=delta43*den0
          weight(3)=num1/den1
! weight(4)          
          weight(4)=wt-weight(1)-weight(2)-weight(3)

        case default

          print*, 'case default in ksurf.f90'

      end select

      end subroutine ksurf

!EOC        
