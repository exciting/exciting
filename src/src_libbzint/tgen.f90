!BOP
!
! !ROUTINE: tgen
!
! !INTERFACE:
      subroutine tgen(ntet,outet,wtet)
!
! !DESCRIPTION:
! The tetrahedra are tried to be related by the symmetry. The tetrahedron is defined
! by the irreducible k-point on the vertices. If the vertices of two tetrahedra goes
! into the same irreducible point set, the weight of this tetrahedron will be added 
! once to reflect this relation. After this, we can just do the integration over the 
! irreducible tetrahedron and multiply it by a weight before summation. 

! 
! !USES:
      
      use order, only: sort
      use kgen_internals
      use tetra_internal, only: redtet
      
      implicit none
      
! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: outet(4,*) ! Id. number of the vertices of tetrahedron i
      integer(4), intent(out) :: wtet(*)    ! Weight of each tetrahedron
      integer(4), intent(out) :: ntet       ! Number of tetrahedra
      
! !LOCAL VARIABLES:      
 
      integer(4) :: i1 
      integer(4) :: i2 
      integer(4) :: i3 
      integer(4) :: i 
      integer(4) :: itet
      integer(4) :: j 
      integer(4) :: t 
      integer(4) :: idx
      integer(4) :: idx2
      integer(4) :: orig(3) 
      integer(4) :: cornid
      integer(4) :: corn(3) 
      integer(4), dimension(4) :: intet,inx
      integer(4), dimension(3,4,6) :: tet
      integer(4), external :: idkp
      logical    :: notfd

! !SYSTEM ROUTINES:
      intrinsic mod
      
      external tetinit
 
!EOP
!BOC
!
!     Initialize the tetrahedra, and the maximum number of tetrahedra
!
      call tetinit(tet)
      ntet=6*div(1)*div(2)*div(3)
      allocate(redtet(ntet))
!
!     Calculate the tetrahedra volume
!
      vt=1.d0/dble(ntet)

      outet(1:4,1:ntet)=0
      idx=0
      idx2=0
      do i1=0,div(1)-1
        orig(1)=i1
        do i2=0,div(2)-1
          orig(2)=i2
          do i3=0,div(3)-1
            orig(3)=i3
            idx2=idx2+1
            do t=1,6
              do i=1,4
                do j=1,3
                  corn(j)=mod(orig(j)+tet(j,i,t),div(j))
                enddo
                cornid=idkp(corn)
                intet(i)=ikpid(redkp(cornid))
              enddo
              call sort(4,intet,inx)
              
              itet=0
              notfd=.true.
              do while ((itet.lt.idx).and.notfd) 
                itet=itet+1
                if((outet(1,itet).eq.intet(1)).and.   &
     &             (outet(2,itet).eq.intet(2)).and.   &
     &             (outet(3,itet).eq.intet(3)).and.   &
     &             (outet(4,itet).eq.intet(4)))then
                  wtet(itet)=wtet(itet)+1
                  redtet(6*(idx2-1)+t)=itet
                  notfd=.false.
                endif
              enddo
              if(notfd)then
                idx=idx+1
                do i=1,4
                  outet(i,idx)=intet(i)
                enddo
                wtet(idx)=1
                redtet(6*(idx2-1)+t)=idx
              endif
            enddo
          enddo
        enddo
      enddo
      ntet=idx

      end subroutine tgen
!EOC
