
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

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
!<sag>
      use control, only: tetraifc, tetradbglv
!</sag>
      implicit none
      
! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: outet(4,*) ! Id. number of the vertices
!                                             of tetrahedron i
     
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
      integer(4) :: index
      integer(4) :: index2
      integer(4) :: orig(3) 
      integer(4) :: cornid
      integer(4) :: corn(3) 
      integer(4), dimension(4) :: intet,inx
      integer(4), dimension(3,4,6) :: tet
      integer(4), external :: idkp
      logical    :: notfd

!!!!!!!!!!!! SAG
integer :: nkptnr

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
      index=0
      index2=0
!<sag>
nkptnr=div(1)*div(2)*div(3)
      if (trim(tetraifc)=='wien2k') then
      ! original code

      do i1=0,div(1)-1
        orig(1)=i1
        do i2=0,div(2)-1
          orig(2)=i2
          do i3=0,div(3)-1
            orig(3)=i3
            index2=index2+1
            do t=1,6
              do i=1,4
                do j=1,3
                  corn(j)=mod(orig(j)+tet(j,i,t),div(j))
                enddo
                cornid=idkp(corn)
                intet(i)=ikpid(redkp(cornid))
!                intet(i)=cornid
              enddo
              call sort(4,intet,inx)
              
              itet=0
              notfd=.true.
              do while ((itet.lt.index).and.notfd) 
                itet=itet+1
                if((outet(1,itet).eq.intet(1)).and.   &
     &             (outet(2,itet).eq.intet(2)).and.   &
     &             (outet(3,itet).eq.intet(3)).and.   &
     &             (outet(4,itet).eq.intet(4)))then
                  wtet(itet)=wtet(itet)+1
                  redtet(6*(index2-1)+t)=itet
                  notfd=.false.
                endif
              enddo
              if(notfd)then
                index=index+1
                do i=1,4
                  outet(i,index)=intet(i)
                enddo
                wtet(index)=1
                redtet(6*(index2-1)+t)=index
              endif
            enddo
          enddo
        enddo
      enddo

      ! end original code
      else if (trim(tetraifc)=='exciting') then
      ! new code

      do i3=0,div(3)-1
        orig(3)=i3
        do i2=0,div(2)-1
          orig(2)=i2
          do i1=0,div(1)-1
            orig(1)=i1
            index2=index2+1
            do t=1,6
              do i=1,4
                do j=1,3
                  corn(j)=mod(orig(j)+tet(j,i,t),div(j))
                enddo
                cornid=idkp(corn)
                intet(i)=ikpid(redkp(cornid))
!                intet(i)=cornid
              enddo
              call sort(4,intet,inx)
              
              itet=0
              notfd=.true.
              do while ((itet.lt.index).and.notfd) 
                itet=itet+1
                if((outet(1,itet).eq.intet(1)).and.   &
     &             (outet(2,itet).eq.intet(2)).and.   &
     &             (outet(3,itet).eq.intet(3)).and.   &
     &             (outet(4,itet).eq.intet(4)))then
                  wtet(itet)=wtet(itet)+1
                  redtet(6*(index2-1)+t)=itet
                  notfd=.false.
                endif
              enddo
              if(notfd)then
                index=index+1
                do i=1,4
                  outet(i,index)=intet(i)
                enddo
                wtet(index)=1
                redtet(6*(index2-1)+t)=index
              endif
            enddo
          enddo
        enddo
      enddo

      ! end new code
      end if ! if (tetraifc)
!</sag>
      ntet=index


      ! *** DEBUG ***
      if (tetradbglv.gt.1) then
         write(*,*) 'TGEN REPORTS:'
         write(*,*) 'nkptnr:',nkptnr
         write(*,*) 'maxtet:',6*nkptnr
         write(*,*) 'ntet  :',ntet
         do j=1,nkptnr
            do itet=1,6
               i=(j-1)*6+itet
               write(*,'(3i5,3x,i9,3x,4i5)') j,itet,i,redtet(i),&
                    outet(:,redtet(i))
            end do
         end do
      end if      

      end subroutine tgen
!EOC
