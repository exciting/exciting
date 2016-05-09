!BOP
! 
! !ROUTINE: cartezian
!
! !INTERFACE:
      subroutine cartezian(nkp,aaa,br2,klist,idiv)

! !DESCRIPTION:
!
! This subroutine transform the submesh coordinates of the kpoints into
! cartesian coordinates in the reciprocal space
!
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: nkp      ! Number of k-points
      real(8),    intent(in) :: aaa(3)   ! lattice constants
      real(8),    intent(in) :: br2(3,3) ! reciprocal lattice vectors
      
! !INPUT/OUTPUT PARAMETERS:
      integer(4), intent(inout) :: klist(3,*) ! integer coordinates of the kpoints      
      integer(4), intent(inout) :: idiv       ! minimum common divisor for the integer coordinates of the kpoints
       
! !LOCAL VARIABLES:

      integer(4) :: kpi
      integer(4) :: i
      integer(4) :: j
      real(8), dimension(3) :: ak
      
! !DEFINED PARAMETERS:

      real(8), parameter :: pi = 3.14159265358979323846      
      
!EOP      
!
!BOC
      
      do kpi=1,nkp
        do i=1,3
          ak(i)=0.0d0
          do j=1,3
            ak(i)=ak(i)+br2(i,j)*dble(klist(j,kpi))
          enddo
        enddo
        do i=1,3
          klist(i,kpi)=nint(ak(i)/2./pi*aaa(i))
        enddo
      enddo
      call divisi(nkp,idiv,klist)

      end subroutine cartezian

!EOC
