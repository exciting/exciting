!BOP
!
! !ROUTINE: writeqgen
!
! !INTERFACE: 
      subroutine writeqgen

! !DESCRIPTION:
!
! This subroutine writes the identification number of the nodal points of
! each tetrahedron 
!
! !USES:
      
      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: i,j,k,itet,kmax
!EOP
!BOC 

      write(99,*) 
      
      call boxmsg(99,'=',"Irreducible part:")

      write(99,*) "#writeqgen: Nodal points of tetrahedron"
  
      write(99,101) nkpt,ntet,tvol
      do itet=1,ntet
        write(99,*)itet,(tnodes(i,itet),i=1,4),wtet(itet)
      enddo 

      call boxmsg(99,'=',"Non reducible part:")
      
      write(99,*) "#writeqgen: Nodal points of tetrahedron"
  
      write(99,101) nkptnr,ntetnr,tvol
      do itet=1,ntetnr
        write(99,*)itet,(tnodesnr(i,itet),i=1,4),wtetnr(itet),(linkq(itet,j),j=1,nkptnr)
      enddo 
!
!     Write the k-dependent q and k' weights
!
      write(99,*)' k,k-q links'
      do i=1,nkptnr,8
        kmax=7
        if((nkptnr-i).lt.kmax)kmax=nkptnr-i
        write(99,103)(i+k,k=0,kmax)
        do j=1,nkptnr
          write(99,102)j,(kqid(j,i+k),k=0,kmax)
        enddo
      enddo

  100 format(i6,4i4,4i6)
  101 format(2i6,e16.8)
  102 format(i4,'   |',8i4)      
  103 format(' ik/iq |',8i4)      
      
      end subroutine writeqgen
!EOC
