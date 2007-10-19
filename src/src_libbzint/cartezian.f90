!BOP
! 
! !ROUTINE: cartezian
!
! !INTERFACE:
      subroutine cartezian(nkp,div,aaa,rbas,klist,idiv)

! !DESCRIPTION:
!
! This subroutine transform the submesh coordinates of the kpoints into
! cartesian coordinates in the reciprocal space
!
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: nkp      ! Number of k-points
      
      integer(4), intent(in) :: div(3)   ! number of divisions of the 
!                                          submesh in each direction

      real(8),    intent(in) :: aaa(3)   ! lattice constants
      
      real(8),    intent(in) :: rbas(3,3)  ! Basis vectors of the direct
!                                          lattice
      
! !INPUT/OUTPUT PARAMETERS:

      integer(4), intent(inout) :: klist(3,*)! integer coordinates of 
!                                              the kpoints      

      integer(4), intent(inout) :: idiv     ! minimum common divisor 
!                                             for the integer 
!                                             coordinates of the
!                                             kpoints      
       
! !LOCAL VARIABLES:

      integer(4) :: kpi
      
      integer(4) :: i
      
      integer(4) :: j
      
      real(8), dimension(3) :: ak
      
      real(8), dimension(3,3) :: gbas ! Basis vectors of the reciprocal
!                                       lattice
      
      
! !DEFINED PARAMETERS:

      real(8), parameter :: pi = 3.14159265358979323846      
      
      external gbass

!EOP      
!
!BOC

      call gbass(rbas,gbas)
      
      do kpi=1,nkp
        do i=1,3
          ak(i)=0.0d0
          do j=1,3
            ak(i)=ak(i)+gbas(j,i)*dble(klist(j,kpi))
          enddo
        enddo
        do i=1,3
          klist(i,kpi)=nint(ak(i)/2./pi*aaa(i))
        enddo
      enddo
      call divisi(nkp,idiv,klist)

      end subroutine cartezian

!EOC
