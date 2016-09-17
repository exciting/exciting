!BOP
! 
! !ROUTINE: cart2int
!
! !INTERFACE:
      subroutine cart2int(nkp,rbas,alat,kvecs)

! !DESCRIPTION:
!
! This subroutine transform the submesh coordinates of the kpoints in
! cartesian coordinates to that of internal coordinates in the reciprocal space
!
! !INPUT PARAMETERS:
      implicit none
      integer(4), intent(in) :: nkp
      real(8),    intent(in) :: rbas(3,3)  ! Basis vectors of the direct lattice in the row-wise form 
      real(8),    intent(in) :: alat(3)  
      real(8),    intent(inout)::kvecs(3,nkp)

! !LOCAL VARIABLES:
      integer(4) :: ik,i,j
      real(8), dimension(3) :: ak
!EOP      
!
!BOC
      do ik=1,nkp
        do i=1,3
          ak(i)=0.0d0
          do j=1,3
            ak(i)=ak(i)+rbas(j,i)*kvecs(j,ik)/alat(j)
          enddo
        enddo
        kvecs(1:3,ik)=ak
      enddo
      end subroutine cart2int
!EOC
