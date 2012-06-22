!BOP
!
! !ROUTINE: calcbarcsing
!
! !INTERFACE:
      subroutine calcbarcsing
      
! !DESCRIPTION:
!
! This subroutine calculates the singular part of the bare coulomb 
! potential matrix by expanding it in plane waves and then transforming 
! to the mixed basis.
!
! !USES:
      
      use modmain
      use modgw
!
! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: i, j
      real(8) :: sqfpi
      real(8) :: fpi
      real(8) :: tstart, tend
    
!EOP
!BOC
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

      if(allocated(barcs))deallocate(barcs)
      allocate(barcs(matsiz,matsiz))
      if(allocated(sqbarcs))deallocate(sqbarcs)
      allocate(sqbarcs(matsiz,matsiz))
      
      fpi=4.0d0*pi
      sqfpi=sqrt(fpi)
      do i=1,matsiz
        do j=1,matsiz
          barcs(i,j)=fpi*wi0(i)*conjg(wi0(j))
          sqbarcs(i,j)=sqfpi*wi0(i)*conjg(wi0(j))
        enddo
      enddo
      
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'CALCBARCSING')
      
      return
      end subroutine calcbarcsing
!EOC            
     
               
