!BOP
!
! !ROUTINE: screencoul
!
! !INTERFACE:
      subroutine screencoul
      
! !DESCRIPTION:
!
!
!This subroutine calculates the correlation part of the screened coulomb matrix as:
!
!\begin{equation}
!W_{ij}(\vec{q})=\sum_l{\left(\epsilon^{-1}_{il}-\delta_{il}\right)v_{lj}}      
!\end{equation}
!

! !USES:

      use modmain   
      use modgw

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iom
      
      real(8)    :: tstart, tend
      
      complex(8), allocatable :: tmat1(:,:), tmat2(:,:)

! !REVISION HISTORY:
!
! Created 18. april 2005 by RGA
! Revisited July 2011 by DIN
!
!EOP
!BOC
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
      
      allocate(tmat1(matsiz,matsiz))
      allocate(tmat2(matsiz,matsiz))
      allocate(wxc(matsiz,matsiz,nomeg))
      
      do iom=1,nomeg
        tmat1(1:matsiz,1:matsiz)=inveps(1:matsiz,1:matsiz,iom)
        call zhemm('l','u',matsiz,matsiz,zone,tmat1,matsiz,sqbarc, &
       &            matsiz,zzero,tmat2,matsiz)   
        call zhemm('l','u',matsiz,matsiz,zone,sqbarc,matsiz,tmat2, &
       &            matsiz,zzero,tmat1,matsiz)   
        wxc(1:matsiz,1:matsiz,iom)=tmat1(1:matsiz,1:matsiz)
      enddo
      
      deallocate(tmat1)
      deallocate(tmat2)
      deallocate(inveps)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'SCREENCOUL')
      
      end subroutine screencoul
!EOC        
        
