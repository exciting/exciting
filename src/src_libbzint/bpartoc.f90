!BOP
!
! !ROUTINE: bpartoc
!
! !INTERFACE: 
      subroutine bpartoc(eo,eu,efer,weight)
      
! !DESCRIPTION:
!
!This subroutine calculates the contribution of one tetrahedron to the
!convolution weights in the case that both tetrahedra (the one at $\vec{k}$
!and the linked one at $\vec{k}-\vec{q}$) are partially occupied. For the case
!of nnod=7,8, we further use a big region minus a small region to deal with 
!systematic error. This is for the bulk integration case. 
!

! !USES:

      use polyhedron
 
! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in) :: eo(4) ! the band energies at k

      real(8), intent(in) :: eu(4) ! the band energies at k-q

      real(8), intent(in) :: efer   ! the fermi energy

! !OUTPUT PARAMETERS:

      real(8), intent(out) :: weight(4) ! the weight at each corner
       
! !LOCAL VARIABLES:

      integer(4) :: i,inod,ssg
      real(8), dimension(4) :: w1,w2
      real(8), dimension(6,3)  :: nodtemp6
      real(8), dimension(5,3)  :: nodtemp5
      real(8), dimension(4,3)  :: nodtemp4
      external setnodes
      external generictetra
      external genericfunf
      external genericprism

! !REVISION HISTORY:
! 
! Created 22nd. April 2004 by RGA, last rivised 14th. Dec 2004 by XZL
!
!EOP
!BOC
      e(1:4)=eo(1:4)
      f(1:4)=eu(1:4)
      ef=efer   
      ndabc=1
!----------------------------------------------------------------
!    ndabc=1 means the region to be integrated
!----------------------------------------------------------------
      call setnodes
      ssg=0
      select case(nnod)
      case(0)
        weight(1:4)=0.0d0
      case(1)
        weight(1:4)=0.0d0
      case(2)
        weight(1:4)=0.0d0
      case(3)
        weight(1:4)=0.0d0
      case(4)            
        nodtemp4(1:4,1:3)=intnodes(1:4,1:3)
        call generictetra(nodtemp4,weight)
      case(5) 
        do inod=1,5
          nodtemp5(index(inod,1),1:3)=intnodes(inod,1:3)
        enddo
        call genericfunf(nodtemp5,weight)
        deallocate(index)

      case(6)
        do inod=1,6
          nodtemp6(index(inod,1),1:3)=intnodes(inod,1:3)
        enddo   
        call genericprism(nodtemp6,weight)
        deallocate(index)
      case(7:8)
         ndabc=2
!-----------------------------------------------------------------
!  when nnod>6 for the region to be integrated, we use the region 
!  of e1<ef minused by that of e1<ef and e2<ef to cancel part of 
!  the systematic error of linearizatio method
!  ndabc=2 and then call setnodes is to get the e1<ef region
!-----------------------------------------------------------------
         deallocate(index)
         call setnodes       
         select case(nnod)
         case(0)
           weight(1:4)=0.0d0
         case(1)
           weight(1:4)=0.0d0
         case(2)
           weight(1:4)=0.0d0
         case(3)
           weight(1:4)=0.0d0
         case(4)
            nodtemp4(1:4,1:3)=intnodes(1:4,1:3)
            call generictetra(nodtemp4,weight)
         case(5)
           do inod=1,5
             nodtemp5(index(inod,1),1:3)=intnodes(inod,1:3)
           enddo
           call genericfunf(nodtemp5,weight)
           deallocate(index)
         case(6)
           do inod=1,6
             nodtemp6(index(inod,1),1:3)=intnodes(inod,1:3)
           enddo
           call genericprism(nodtemp6,weight)
           deallocate(index)

         case default
           print*, nnod
           stop 'error in bpartoc.f90(0)'

         end select
         ndabc=3
!----------------------------------------------------------------------
!  ndabc=3 and then setnodes is to get the e1<ef and e2<ef region
!----------------------------------------------------------------------
         call setnodes
         select case(nnod)
         case(0)
           w1(1:4)=0.0d0
         case(1)
           w1(1:4)=0.0d0
         case(2)
           w1(1:4)=0.0d0
         case(3)
           w1(1:4)=0.0d0
         case(4)
           nodtemp4(1:4,1:3)=intnodes(1:4,1:3)
           call generictetra(nodtemp4,w1)
         case(5)
           do inod=1,5
             nodtemp5(index(inod,1),1:3)=intnodes(inod,1:3)
           enddo
           call genericfunf(nodtemp5,w1)
           deallocate(index)
         case(6)
           do inod=1,6
             nodtemp6(index(inod,1),1:3)=intnodes(inod,1:3)
           enddo
           call genericprism(nodtemp6,w1)
           deallocate(index)
         case(7)
           do inod=1,7
             if(index(inod,1).ne.0) then
               nodtemp5(index(inod,1),1:3)=intnodes(inod,1:3)
             else
              continue
             endif
           enddo
           call genericfunf(nodtemp5,w1)
           do inod=1,7
             if(index(inod,2).ne.0) then
               nodtemp5(index(inod,2),1:3)=intnodes(inod,1:3)
             else
               continue
             endif
           enddo
           call genericfunf(nodtemp5,w2)
           do i=1,4
             w1(i)=w1(i)+w2(i)
           enddo
           deallocate(index)
            
         case(8)
           do inod=1,8
             if(index(inod,1).ne.0) then
               nodtemp6(index(inod,1),1:3)=intnodes(inod,1:3)
             else
               continue
             endif
           enddo
           call genericprism(nodtemp6,w1)
           do inod=1,8
             if(index(inod,2).ne.0) then
               nodtemp6(index(inod,2),1:3)=intnodes(inod,1:3)
             else
               continue
             endif
           enddo
           call genericprism(nodtemp6,w2)
           do i=1,4
             w1(i)=w1(i)+w2(i)
           enddo
           deallocate(index)
         case default
            print*, nnod
            stop 'error in bpartoc.f90(1)'
         end select
         do i=1,4
            weight(i)=weight(i)-w1(i)
         enddo
      case default
         print*, nnod
         stop 'error in bpartoc.f90(2)'
      end select


      end subroutine bpartoc

!EOC        
