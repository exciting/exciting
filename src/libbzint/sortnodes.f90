
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: sortnodes
!
! !INTERFACE:
      subroutine sortnodes
      
! !DESCRIPTION:
! This subroutine sort the nodes of each region so that it can be devided
! into tetrahedron in order. For the case of nnod=4, it is already done.
! For nnod=5, we devide it into two tetrahedrons with index(inod,1:2) indicates
! the order of the node in the tetrahedron. For example, if index(inod,2)=1,
! we should take this node as the first node in the second tetrahedron. 
! For nnod=6, we can devide it into three tetrahedron, the nodes should be ordered
! in a way that 1,2,3,4 form a tetrahedron, 2,,3,4,5 form a tetrahedron and 
! 3,4,5,6 form the other one. These three together make up the prism. 
! For nnod=7, we can devide it into two penta which are ordered in the way as
! nnod=5. For nnod=8, we can devide it into two prism which are ordered in the
! way as nnod=6 case.  
!

! !USES:
      
      use polyhedron
      
! !LOCAL VARIABLES:

      implicit none

      integer(1) :: ibit, inod,  ip, iind
      integer(4) :: is, it, i, iq1
      integer(1) :: ibt(3) 
      integer(1) :: it1, it2, is1
      integer(1), dimension(6) :: triangle
      integer(1), dimension(6) :: square
      integer(1), dimension(2) :: penta
      integer(4), dimension(0:5) :: spl
      integer(1) ::  nt, signbit, signbita,signbitb,signnod,signnodb
      integer(1), dimension(6) :: sigtri
      integer(1), parameter :: one=1
      integer(1), parameter :: four=4
      integer(1), parameter :: five=5

! !INTRINSIC ROUTINES:
      
      intrinsic btest  

! !REVISION HISTORY:
!
! Created 23rd. April 2004 by RGA, last revised by XZL on Dec 14th, 2004
!
!EOP
!BOC
      spl(0:5)=0
      it=0
      is=0
      ip=0
      select case(nnod)

      case(5)     
        allocate(index(5,1))
!--------------------------------------------------------------------
!       finding how many triangles and squares to form the region
!                               begin
!--------------------------------------------------------------------

        do ibit=0,5     ! here we use the bit to represent the six planes
          do inod=1,5   ! inod means each node
            nt=ntype(inod)
            if(btest(nt,ibit))spl(ibit)=spl(ibit)+1
          enddo
          select case(spl(ibit))
          case(3)      ! for one plane, if there are 3 nodes on it.
            it=it+1
            triangle(it)=ibit
          case(4)      ! for one plane, if there are 4 nodes on it.
            is=is+1
            square(is)=ibit
            signbit=ibit
          end select
        enddo

!--------------------------------------------------------------------
!       finding how many triangles and squares to form the region
!                               end
!--------------------------------------------------------------------

        index(1:5,1)=0
        if(is.ne.0) then     ! the case when there is one node
!                                     out of one surface while the other
!                                     four in.  
!--------------------------------------------------------------------
!  When the region is formed by 4 triangles and one square, begin
!--------------------------------------------------------------------
         is1=1
         do inod=1,5
           if(btest(ntype(inod),triangle(1))) then
             if(.not.(btest(ntype(inod),signbit))) then
               index(inod,1)=2
             else
               index(inod,1)=is1
               is1=is1+2
               signnod=inod
             endif
           else
             continue
           endif
         enddo
         
         do i=2,4
           if(btest(ntype(signnod),triangle(i)))signbit=triangle(i)
         enddo

         do inod=1,5
           if(index(inod,1).eq.0) then
             if(btest(ntype(inod),signbit)) then
               index(inod,1)=5
             else
               index(inod,1)=4
             endif
           else
             continue
           endif
         enddo
!-------------------------------------------------------------------------
!   when the region is formed by 4 triangles and one square, end
!--------------------------------------------------------------------------

        else                  ! for the case when there are 3 nodes in 
!                               each surface.
!-------------------------------------------------------------------------
!   when the region is formed by 6 triangles, begin
!--------------------------------------------------------------------------

         do inod=1,5
           signbit=0
           do ibit=0,5
             if(btest(ntype(inod),ibit))signbit=signbit+1
           enddo
           it1=1
           it2=1
           if(signbit.eq.3) then
             index(inod,1)=1+(it1-1)*4
             it1=it1+1
           else
             index(inod,1)=1+it2
             it2=it2+1
           endif
         enddo
!-------------------------------------------------------------------------
!   when the region is formed by 6 triangles, end
!--------------------------------------------------------------------------


        endif             

 
      case(6)
        allocate(index(6,1))
        index(1:6,1)=0
!--------------------------------------------------------------------
!       finding how many triangles and squares to form the region
!                               begin
!--------------------------------------------------------------------
       
        do ibit=0,5
          do inod=1,6
            nt=ntype(inod)
            if(btest(nt,ibit))spl(ibit)=spl(ibit)+1
          enddo
          select case(spl(ibit))
          case(3)
            it=it+1
            triangle(it)=ibit
          
          case(4)
            is=is+1
            square(is)=ibit
        
          end select  
        enddo
!--------------------------------------------------------------------
!       finding how many triangles and squares to form the region
!                               end
!--------------------------------------------------------------------
        select case(is)

        case(2)
!-------------------------------------------------------------------------
!   when the region is formed by 4 triangles and 2 squares, begin
!-------------------------------------------------------------------------

          it=1
          do inod=1,6
            if(btest(ntype(inod),square(1)).and.                    &
     &         btest(ntype(inod),square(2))) then
               index(inod,1)=2+(it-1)*3
               if(it.eq.1) signnod=inod
               if(it.eq.2) signnodb=inod               
               it=it+1
            else
               continue
            endif
          enddo
          sigtri(1:6)=0
          do inod=1,6
            do it=1,4
              if(btest(ntype(inod),triangle(it))) then
                 sigtri(inod)=sigtri(inod)+1
              else
                 continue
              endif
            enddo
          enddo
          do it=1,4
            if(btest(ntype(signnod),triangle(it)))signbita=it
            if(btest(ntype(signnodb),triangle(it)))signbitb=it
          enddo
          do inod=1,6
             if(index(inod,1).eq.0) then
                if((sigtri(inod).eq.3).and.btest(ntype(inod),triangle(signbita)))index(inod,1)=4
                if((sigtri(inod).eq.2).and.btest(ntype(inod),triangle(signbita)))index(inod,1)=1
                if((sigtri(inod).eq.3).and.btest(ntype(inod),triangle(signbitb)))index(inod,1)=3
                if((sigtri(inod).eq.2).and.btest(ntype(inod),triangle(signbitb)))index(inod,1)=6
             else
                continue
             endif
          enddo
!-------------------------------------------------------------------------
!   when the region is formed by 4 triangles and 2 squares, end
!-------------------------------------------------------------------------

        case(3)
          do inod=1,6
           do it=1,2
            if(btest(ntype(inod),triangle(it)))then
                if(btest(ntype(inod),square(1)).and.                 &
     &             btest(ntype(inod),square(2))) then
                    index(inod,1)=1+3*(it-1)
                else if(btest(ntype(inod),square(1))) then
                    index(inod,1)=2+3*(it-1)
                else if(btest(ntype(inod),square(2))) then
                    index(inod,1)=3+3*(it-1)
                else
                   stop 'error in sortnodes, nnod=6'
                endif
            else
             continue
            endif
           enddo
          enddo
        end select
      case(7) 
        allocate(index(7,2))
!-----------------------------------------------------------------------
!   Finding how the region is formed, begin
!-----------------------------------------------------------------------
        index(1:7,1:2)=0
        do ibit=0,5
          do inod=1,7
             if(btest(ntype(inod),ibit))spl(ibit)=spl(ibit)+1
          enddo
          select case(spl(ibit))
          case(3)
            it=it+1
            triangle(it)=ibit
          case(4)
            is=is+1
            square(is)=ibit
          case(5)
            ip=ip+1
            penta(ip)=ibit
          end select
        enddo
        write(90,'(a7,3i4)')'istp = ',it,is,ip
!-----------------------------------------------------------------------
!    Finding how the region is formed, end
!-----------------------------------------------------------------------

        select case(it)
        case(2)       ! when it is formed by 2 triangles and 4 squares
!-----------------------------------------------------------------------
! When the region is formed by 2 triangles and 4 squares, begin
!-----------------------------------------------------------------------
          select case(is)
          case(4)    ! 2 triangles and 4 squares
            do inod=1,7
             if(btest(ntype(inod),triangle(1)).and.  &
     &          btest(ntype(inod),triangle(2))) then
                index(inod,1:2)=2
                signnod=inod
             else
                continue
             endif
            enddo
     
            is1=0
            do ibit=0,5
              if(btest(ntype(signnod),ibit).and.  &
     &           (spl(ibit).eq.4).and.(is1.eq.0)) then
                signbita=ibit
                is1=is1+1
              else if(btest(ntype(signnod),ibit).and. &
     &           (spl(ibit).eq.4)) then
                signbitb=ibit
              else
                continue
              endif
            enddo
            do inod=1,7
             if(inod.ne.signnod) then 
               if(btest(ntype(inod),triangle(1)).and.  &
     &            btest(ntype(inod),signbita)) then
                  index(inod,1)=1
               else if(btest(ntype(inod),triangle(1)).and.  &
     &            btest(ntype(inod),signbitb)) then
                  index(inod,1)=3
               else if(btest(ntype(inod),triangle(2)).and.  &
     &            btest(ntype(inod),signbita)) then
                  index(inod,2)=1
               else if(btest(ntype(inod),triangle(2)).and.  &
     &            btest(ntype(inod),signbitb)) then
                  index(inod,2)=3
               else if(btest(ntype(inod),signbita)) then
                  index(inod,1:2)=4
               else
                  index(inod,1:2)=5
               endif
             else
               continue
             endif
            enddo
          case default
            print*, it, is, ip
            stop 'error in sortnodes.f90, inod=7 and it=2'
          end select
!-----------------------------------------------------------------------
! When the region is formed by 2 triangles and 4 squares, end
!-----------------------------------------------------------------------

        case(3)  ! when it is made up of 3 triangles, 2 squares and 1 penta
!--------------------------------------------------------------------------
! When the region is formed by 3 triangles and 2 squares and 1 penta, begin
!--------------------------------------------------------------------------
          index(1:7,1:2)=0
          do inod=1,7
            ibit=0
            do it1=1,3
              if(btest(ntype(inod),triangle(it1)))ibit=ibit+1
            enddo
            if(btest(ntype(inod),penta(1))) then
              if(ibit.eq.2)index(inod,1)=1
            else
              if(ibit.eq.2) then
                index(inod,1)=2
                index(inod,2)=1
              else
                index(inod,2)=4
                iind=inod
              endif
            endif
          enddo
          do it2=1,3
            if(btest(ntype(iind),triangle(it2)))signbita=triangle(it2)
          enddo
          do inod=1,7
            if((index(inod,1).eq.0).and.(index(inod,2).eq.0)) then
              if(btest(ntype(inod),signbita).and.  &
     &           btest(ntype(inod),square(1))) then
                 index(inod,2)=5
              else if(btest(ntype(inod),signbita).and.  &
     &           btest(ntype(inod),square(2))) then
                 index(inod,1)=5
                 index(inod,2)=2
              else if(btest(ntype(inod),square(1))) then
                 index(inod,1)=4
                 index(inod,2)=3
              else 
                 index(inod,1)=3
              endif
            else
              continue
            endif
          enddo

        case default
          print*, it, is, ip
          stop 'error in sortnodes.f90, inod=7'
        end select
!--------------------------------------------------------------------------
! When the region is formed by 3 triangles and 2 squares and 1 penta, begin
!--------------------------------------------------------------------------

      case(8)
        allocate(index(8,2))
!--------------------------------------------------------------------------
!  Finding how the region is formed, begin
!--------------------------------------------------------------------------
        do ibit=0,5
          do inod=1,8
             nt=ntype(inod)
             if(btest(nt,ibit))spl(ibit)=spl(ibit)+1
          enddo

          select case(spl(ibit))
          case(3)
            it=it+1
            triangle(it)=ibit
          case(4)
            is=is+1
            square(is)=ibit
          case(5)
            ip=ip+1
            penta(ip)=ibit
          end select
        enddo
!--------------------------------------------------------------------------
! Finding how the region is formed, end
!--------------------------------------------------------------------------    
        index(1:8,1:2)=0
!        write(90,*)'itsp =',it, is,ip
        select case(is)
        case(2)         ! then it is made up of 2 squares, 2 triangles and 2 penta
!--------------------------------------------------------------------------
! When the region is formed by 2 triangles and 2 squares and 2 penta, begin
!--------------------------------------------------------------------------
         do inod=1,8
           do it1=1,2
             if(btest(ntype(inod),triangle(it1))) then
               if(btest(ntype(inod),penta(1)).and. &
     &            btest(ntype(inod),penta(2)))  then
                  index(inod,1)=1+(it1-1)*3
               else if(btest(ntype(inod),penta(1))) then
                  index(inod,1)=2+(it1-1)*3
               else if(btest(ntype(inod),penta(2))) then
                  index(inod,1)=3+(it1-1)*3
               else
                  continue
               endif
             else
               continue
             endif
           enddo
        
           do iq1=1,2
             if(btest(ntype(inod),penta(iq1))) then
               if(btest(ntype(inod),square(1)).and.  &
     &            btest(ntype(inod),square(2))) then
                  index(inod,2)=1+(iq1-1)*3
               else if(btest(ntype(inod),square(1))) then
                  index(inod,2)=2+(iq1-1)*3
               else if(btest(ntype(inod),square(2))) then
                  index(inod,2)=3+(iq1-1)*3
               else
                  continue
               endif
             else
               continue
             endif
           enddo  

         enddo         
!--------------------------------------------------------------------------
! When the region is formed by 2 triangles and 2 squares and 2 penta, end
!--------------------------------------------------------------------------

        case(6)      ! when it is made up of 6 squares

!--------------------------------------------------------------------------
! When the region is formed by 4 squares, begin
!--------------------------------------------------------------------------
         is1=0

! choose the node 1 to be 1 in the first prism and it is not in the 
! second (node 1:0)
         index(1,1)=1
         index(1,2)=0

! set the three planes to which node 1 belongs as 1,2,3 in order of
! appearence         
         do ibit=0,5
           if(btest(ntype(1),ibit)) then
             is1=is1+1
             ibt(is1)=ibit
           endif
         enddo
! the rest of the nodes         
         do inod=2,8
         
           if(btest(ntype(inod),ibt(1)))then
             if(btest(ntype(inod),ibt(2)))then
               index(inod,1:2)=2 ! node 2:2= intersection of plane 1 and 2
             else if(btest(ntype(inod),ibt(3)))then
               index(inod,1:2)=3 ! node 3:3= intersection of plane 1 and 4
             else
               index(inod,1)=0  ! node 0:1= the remaining node of plane 1
               index(inod,2)=1
             endif
           else if(btest(ntype(inod),ibt(2)))then
             if(btest(ntype(inod),ibt(3)))then
               index(inod,1)=4  ! node 4:0= intersection of plane 2 and 3
               index(inod,2)=0
             else
               index(inod,1:2)=5 ! node 5:5= the remaining node of plane 2
             endif
           else
             if(btest(ntype(inod),ibt(3)))then
               index(inod,1:2)=6 ! node 6:6= the remaining node of plane 3
             else
               index(inod,1)=0   ! node 0:4= the node that doesn´t belong
               index(inod,2)=4   ! to any of the planes 1,2 or 3
             endif
           endif
         enddo              
!--------------------------------------------------------------------------
! When the region is formed by 4 squares, end
!--------------------------------------------------------------------------
        case default

          stop 'error in nnod=8 of sortnodes.f90'         

        end select
            
      end select
            
      end subroutine sortnodes
!EOC            
              
            
      
