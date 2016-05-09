!BOP
!
! !ROUTINE: unrepnodes
!
! !INTERFACE:
      subroutine unrepnodes(iwrite)
      
! !DESCRIPTION:
!
! This subroutines select, from the intersections of the planes, only
!those that belong to the surface of the tetrahedron. If one of them is
!repeated, which means that there are more than three planes intersecting,
!then it set the corresponding bits of \verb{ntype} to 1.
!

! !USES:
      use tetra_internal, only: fout
      use polyhedron

! !LOCAL VARIABLES:

      implicit none
      integer(4), intent(in) :: iwrite
      integer(1) :: i,j,k,ib
      integer(1) :: bti,btj
      integer(4) :: inod
      real(8), dimension(3) :: intnodi
      real(8) :: sumnt
      integer(1), parameter :: one=1
      real(8), parameter :: zerotol = 1.0d-8
      
      intrinsic ibits
      intrinsic ibset 

! !REVISION HISTORY:
!
! Created 22nd. April 2004 by RGA
!
!EOP
!BOC
      i = 1
      do while (i.lt.nnod)
        intnodi(1:3)=intnodes(1:3,i)
        j=i+1
        do while (j.le.nnod)
          sumnt=0.0d0
          do k=1,3
            sumnt=sumnt+dabs(intnodi(k)-intnodes(k,j))   ! this is mainly to 
!                                    see whether these two points are too close
          enddo
          if(sumnt.lt.zerotol)then                    ! if they are close enough
!                                   to be viewed as the same point, delete one of
!                                   them and keep the other.
            do ib=0,5
              btj=ibits(ntype(j),ib,one)
              bti=ibits(ntype(i),ib,one)
              if((btj.eq.1).and.(bti.eq.0))ntype(i)=ibset(ntype(i),ib)
            enddo 
            do k=j+1,nnod
              intnodes(1:3,k-1)=intnodes(1:3,k)
              ntype(k-1)=ntype(k)
            enddo
            intnodes(1:3,nnod)=0.0d0
            ntype(nnod)=0
            nnod=nnod-1
        else
           j=j+1
        endif
       enddo
       i=i+1
      enddo
      if(iwrite.eq.1)then
        write(fout,*)'unrepeated nodes (from unrepnodes)'
        do inod=1,nnod
          write(fout,'(i4,3f18.10,i4,1x,b6.6)')inod, intnodes(1:3,inod),ntype(inod),ntype(inod)
        enddo
        write(fout,*)
      endif    

      end subroutine unrepnodes
!EOC      
