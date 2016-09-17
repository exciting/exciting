!BOP
!
! !ROUTINE: tetinit
!
! !INTERFACE:
      subroutine tetinit(tet)
       
! !DESCRIPTION:
!
! This subroutine selects the tetrahedra so that they share the shortest 
! diagonal of the containing cube. The tetrahedra will be divided according
! to this diagonal to lower the possible error of linearization method. 
!
      
! !USES:

      use kgen_internals      
      use tetra_internal, only: mndg

      implicit none

! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: tet(3,4,6)

! !LOCAL VARIABLES:
      integer(4) :: i,j,k
      real(8), dimension(4) :: diag
      real(8) :: a

! !DEFINED PARAMETERS:

      integer p0(8,3),tet0(3,4,6)
      data tet0 / 0,0,0, 0,0,1, 0,1,1, 1,1,1, &
     &            0,0,0, 0,1,1, 0,1,0, 1,1,1, &
     &            0,0,0, 0,1,0, 1,1,0, 1,1,1, &
     &            0,0,0, 1,1,0, 1,0,0, 1,1,1, &
     &            0,0,0, 1,0,0, 1,0,1, 1,1,1, &
     &            0,0,0, 1,0,1, 0,0,1, 1,1,1/
      data p0/ 0,0,0,0,1,1,1,1, &
     &         0,0,1,1,0,0,1,1, &
     &         0,1,0,1,0,1,0,1/

!EOP
!BOC

! calculate main diagonals
      
      do i=1,4
        diag(i)=0.d0
        do j=1,3
          a=0.d0
          do k=1,3
            a=a+gbas(k,j)*(p0(i,k)-p0(9-i,k))/dble(div(k))
          enddo
          diag(i)=diag(i)+a**2
        enddo
      enddo

! find smallest diagonal
      mndg=1
      do i=2,4
        if(diag(i).lt.diag(mndg)) mndg=i
      enddo

! rotate tetraedra
      do i=1,6
      do j=1,4
        if(mndg.eq.1) then
          tet(3,j,i)=  tet0(3,j,i)
          tet(2,j,i)=  tet0(2,j,i)
          tet(1,j,i)=  tet0(1,j,i)
        endif
        if(mndg.eq.2) then
          tet(3,j,i)=1-tet0(2,j,i)
          tet(2,j,i)=  tet0(3,j,i)
          tet(1,j,i)=  tet0(1,j,i)
        endif
        if(mndg.eq.3) then
          tet(3,j,i)=  tet0(2,j,i)
          tet(2,j,i)=1-tet0(3,j,i)
          tet(1,j,i)=  tet0(1,j,i)
        endif
        if(mndg.eq.4) then
          tet(3,j,i)=1-tet0(3,j,i)
          tet(2,j,i)=1-tet0(2,j,i)
          tet(1,j,i)=  tet0(1,j,i)
        endif
      enddo
      enddo
         
      return
      end subroutine tetinit
!EOC
