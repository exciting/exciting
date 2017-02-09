!BOP
!
! !ROUTINE: relnodes
!
! !INTERFACE:
      subroutine relnodes(iwrite,info)
!     
! !DESCRIPTION:
!
! This subroutine select the relevant intnodes, that is, only those defining
! the polyhedron that has to be integrated. The intnodes selected are those for
! which $\varepsilon(\vec{k})\le \varepsilon_F$ and 
! $\varepsilon'(\vec{k})\ge \varepsilon_F$ if ndabc=1. $\varepsilon(\vec{k})
! \le \varepsilon_F$ when ndabc=2.  $\varepsilon(\vec{k})\le \varepsilon_F$ and
! $\varepsilon'(\vec{k})\le \varepsilon_F$ when ndabc=3. So the first region equals
! the second region minused by the third region.

!
! !USES:
      use tetra_internal, only: fout
      use polyhedron
      
! !INPUT PARAMETERS:      

      implicit none
      
      integer(4), intent(in) :: iwrite ! flag for writing extra output

! !OUTPUT PARAMETERS:
      
      integer(4), intent(out) :: info ! flag for correct execution      

! !LOCAL VARIABLES:

      integer(4) :: inod, jnod, maxnod
      
      real(8) :: e1, e2
      real(8), parameter :: zerotol=1.00d-12 
      
      real(8), dimension(3) :: node_coord 
      
      logical :: keep_node, e1_occup, e2_unocc
      
      real(8), external :: tlinap

! !REVISION HISTORY:
!
! Created 22nd. April 2004 by RGA
! Last revised 14th. Dec 2004 by XZL
!
!EOP
!BOC
      info = 0
      inod = 1
      maxnod = nnod
      if(iwrite.eq.1)then
        write(fout,*)'energies at the tetracorn'
        write(fout,'("e_occ   =",4f18.10)')e
        write(fout,'("e_unocc =",4f18.10)')f
      endif  
      do while (inod.le.maxnod)
        node_coord(1:3)=intnodes(1:3,inod)
        e1=tlinap(node_coord,e)-ef
        e2=tlinap(node_coord,f)-ef
        
        select case (ntype(inod))
!
! The corners of the tetrahedron, both conditions have to be satisfied        
!
          case(7:15) 
            e1_occup = (e1.lt.0.0d0).or.(dabs(e1).le.zerotol)
            e2_unocc = (e2.gt.0.0d0).or.(dabs(e2).le.zerotol)
            keep_node = e1_occup .and. e2_unocc
!
! Intersections of plane 5 (e1=ef) with two planes of the tetrahedron
! only e2_unocc has to be satisfied.
!
          case(16:31)
            e2_unocc = (e2.gt.0.0d0).or.(dabs(e2).le.zerotol)
            keep_node = e2_unocc
!
! Intersections of plane 6 (e2=ef) with two planes of the tetrahedron
! only e1_occup has to be satisfied.
!
          case(32:47)
            e1_occup = (e1.lt.0.0d0).or.(dabs(e1).le.zerotol)
            keep_node = e1_occup
!
! Intersections of plane 5 (e1=ef) and 6 (e2=ef) with one or more planes 
! of the tetrahedron: Both conditions are already fullfiled
!
          case(48:63)
            keep_node=.true.
!
! The value of ntype is wrong, write message and send error signal.
!          
          case default
            write(fout,*)'ERROR: ntype(inod) =',ntype(inod),' > 63 !!!'
            info = 1  

        end select         

        if(iwrite.eq.1)write(fout,'(i4,2g18.10,l2)')inod,e1,e2,keep_node
 
        if(keep_node)then  
          inod=inod+1
        else ! remove the node from the list
          do jnod=inod+1,maxnod
            intnodes(1:3,jnod-1)=intnodes(1:3,jnod)
            ntype(jnod-1)=ntype(jnod)
          enddo
          maxnod=maxnod-1
        endif
      enddo ! inod

      if(iwrite.eq.1)call flushbuf(6)
!      
! Set the variables of the removed nodes to zero, to avoid errors.
!
      do inod=maxnod+1,nnod
        intnodes(1:3,inod)=0.0d0
        ntype(inod)=0
      enddo
! 
! Set the number of nodes to the relevant ones
!      
      nnod=maxnod

      if(iwrite.eq.1)then
        write(fout,*)'relevant nodes (from relnodes)'
        do inod=1,nnod
          write(fout,'(i4,3f18.10,i4,1x,b6.6)')inod, intnodes(1:3,inod),ntype(inod),ntype(inod)
        enddo
        write(fout,*)
        call flushbuf(6)
      endif    

      end subroutine relnodes
      
!EOC 
