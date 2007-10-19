!BOP
!
! !ROUTINE: kqpartoc
!
! !INTERFACE: 
      subroutine kqpartoc(eo,eu,efer,weight)
      
! !DESCRIPTION:
!
!This subroutine calculates the contribution of one tetrahedron to the
!convolution weights in the case that both tetrahedra (the one at $\vec{k}$
!and the linked one at &\vec{k}-\vec{q}$) are partially occupied.
!

! !USES:

      use polyhedron
 
! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in) :: eo(4) ! the band energies at k

      real(8), intent(in) :: eu(4) ! the band energies at k-q

      real(8), intent(in) :: efer   ! the fermi energy

!      real(8), intent(in) :: omeg

!      integer(4), intent(in) :: sigfreq
! !OUTPUT PARAMETERS:

      real(8), intent(out) :: weight(4) ! the weight at each corner
       
! !LOCAL VARIABLES:

      integer(4) :: i,inod
      real(8), dimension(4) :: www,wwt
      real(8), dimension(6,3)  :: nodtemp, nodtemp0, nodtemp1
      real(8), dimension(4,3)  :: nodtemp4
      external setnodes
      external generictetra
      external genericprism

! !REVISION HISTORY:
! 
! Created 22nd. April 2004 by RGA, last rivised 1st. September by XZL
!
!EOP
!BOC
      e(1:4)=eo(1:4)
      f(1:4)=eu(1:4)
      ef=efer
     
      call surfnodeskqmesh
      
      call unrepnodes

      call relnodeskq

      call sortnodes
!     print*, nnod
      select case(nnod)
      case(0)
         wwt(1:4)=0.0d0
      case(1)
         wwt(1:4)=0.0d0
      case(2)
         wwt(1:4)=0.0d0
      case(3)
         wwt(1:4)=0.0d0   
      case(4)            
        nodtemp4(1:4,1:3)=intnodes(1:4,1:3)
        call generictetra(nodtemp4,wwt)
!        print*, wwt(:)      
   
      case(5) 
        do inod=1,5
          nodtemp(inod,1:3)=intnodes(index(inod,1),1:3)
!          print*, nodtemp(inod,:)
        enddo
        call genericfunf(nodtemp,wwt)
      
      case(6)
        do inod=1,6
          nodtemp(inod,1:3)=intnodes(index(inod,1),1:3)
        enddo   
        call genericprism(nodtemp,wwt)
!        print*,wwt(:)
      case default
          stop 'error in kpartoc.f90'
 
      end select

      weight(1:4)=wwt(1:4)
      deallocate(index)
      deallocate(index0)
      deallocate(index1)


      end subroutine kqpartoc

!EOC        
