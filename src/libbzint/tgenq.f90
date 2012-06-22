
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: tgenq
!
! !INTERFACE:
      subroutine tgenq(q,maxtet,outet,sib)
!
! !DESCRIPTION:
! Given a q vector in the submesh coordinates, this subroutine will give you the
! data about how one tetrahedron is related with another one by this q vector.
! outet has nothing to do with the q vector, it just gives the number k-points on 
! the vertices of each tetrahedron. While sib(:) tells us how one tetrahedron is 
! related to another by the q-vector. If sib(1)=10, it means the first tetrahedron
! is related to the tenth tetrahedron by this q-vector. The tetrahedron is not 
! reduced by symmetry in this subroutine. 

!     
! !USES:
      
      use kgen_internals
!<sag>
      use control, only: tetraifc,kplusq
!</sag>

      implicit none      
      
! !INPUT PARAMETERS:
 
      integer(4), intent(in) :: q(3)        ! submesh coordinates of q
     
! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: maxtet    ! equals 6*ndiv(1)*ndiv(2)*ndiv(3). It tells
!                                        us how many tetrahedra we get without symmetry
      integer(4), intent(out) :: outet(4,*)  ! The id number of the 4 vertices for one
!                                              tetrahedron 
      
      integer(4), intent(out) :: sib(*)    ! just linkt(:) in the calling program, it 
!                                       tells us how one tetrahedron is related to another
!                                       with this q vector
      
! !LOCAL VARIABLES:      
 
      integer(4) :: i1 
      integer(4) :: i2 
      integer(4) :: i3 
      integer(4) :: i 
      integer(4) :: j 
      integer(4) :: t 
      integer(4) :: index
      integer(4) :: orig(3) 
      integer(4) :: cornid
      integer(4) :: corn(3) 
      integer(4) :: orig2(3)
      integer(4), dimension(3,4,6) :: tet
      integer(4), external :: idkp

! !SYSTEM ROUTINES:
      intrinsic mod
      
      external tetinit

!EOP
!BOC
      call tetinit(tet)
      maxtet=6*div(1)*div(2)*div(3)
      vt=1.d0/dble(maxtet)
      index=0
!<sag>
      if (trim(tetraifc)=='wien2k') then

      ! original code
      do i1=0,div(1)-1
        orig(1)=i1
        do i2=0,div(2)-1
          orig(2)=i2
          do i3=0,div(3)-1
            orig(3)=i3
            index=index+1
            do t=1,6
              do i=1,4
                do j=1,3
                  corn(j)=mod(orig(j)+tet(j,i,t),div(j))
                enddo
                cornid=idkp(corn)
                outet(i,6*(index-1)+t)=cornid
              enddo
                do j=1,3
                  if (kplusq) then
                     ! k+q
                     orig2(j)=mod(orig(j)+q(j)+(1-isign(1,orig(j)+q(j)))/2*div(j),div(j)) !SAG
                  else
                     ! k-q
                     orig2(j)=mod(orig(j)-q(j)+(1-isign(1,orig(j)-q(j)))/2*div(j),div(j))
                  end if
                enddo
              sib(6*(index-1)+t)=6*(idkp(orig2)-1)+t
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
            index=index+1
            do t=1,6
              do i=1,4
                do j=1,3
                  corn(j)=mod(orig(j)+tet(j,i,t),div(j))
                enddo
                cornid=idkp(corn)
                outet(i,6*(index-1)+t)=cornid
              enddo
                do j=1,3
                  if (kplusq) then
                     ! k+q
                     orig2(j)=mod(orig(j)+q(j)+(1-isign(1,orig(j)+q(j)))/2*div(j),div(j)) !SAG
                  else
                     ! k-q
                     orig2(j)=mod(orig(j)-q(j)+(1-isign(1,orig(j)-q(j)))/2*div(j),div(j))
                  end if
                enddo
              sib(6*(index-1)+t)=6*(idkp(orig2)-1)+t
            enddo
          enddo
        enddo
      enddo
      ! end new code
      end if ! if (tetraifc)
!</sag>

      end subroutine tgenq
!EOC
