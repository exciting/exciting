!BOP
!
! !MODULE: order
! \begin{verbatim}
      module order
      
!\end{verbatim}
!      
! !INTERFACE:
        interface sort
          module procedure sorti
          module procedure sortr
        end interface
      contains
!EOP      
!
!BOP
! !IROUTINE: sorti
!
! !INTERFACE:
      subroutine sorti(n,v,index)
      
! !DESCRIPTION:
!
! This subroutine sorts the values of the vector v (integer) in increasing
!order, optionally, it also returns an integer vector that gives the old
!indexes in vector v according to the new ones
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: n ! dimension of the vector v
      
! !INPUT/OUTPUT PARAMETERS:

      integer(4), intent(inout) :: v(n)
      
      integer(4), optional, intent(inout) :: index(n)

! !LOCAL VARIABLES:

      integer(4) :: i,vtmp,j
      integer    :: itmp
      
! !SYSTEM ROUTINES:
      

! !REVISION HISTORY:
!
! Created 26th. Feb. 2004
!
!EOP
!BOC
      if(present(index))then
        do i=1,n
          index(i)=i
        enddo
      endif  
      do i=1,n-1
        vtmp=v(i)
        itmp=i
        do j=i+1,n
          if(v(j).lt.vtmp)then
            vtmp=v(j)
            itmp=j
          endif  
        enddo
        if(itmp.ne.i)then
          vtmp=v(i)
          v(i)=v(itmp)
          v(itmp)=vtmp
          if(present(index))then
            vtmp=index(i)
            index(i)=index(itmp)
            index(itmp)=vtmp
          endif  
        endif
      enddo
      end subroutine sorti
!EOC                
!BOP
! !IROUTINE: sortr
!
! !INTERFACE:
      subroutine sortr(n,v,index)
      
! !DESCRIPTION:
!
! This subroutine sorts the values of the vector v (real) in increasing 
!order, optionally, it also returns an integer vector that gives the old
!indexes in vector v according to the new ones
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: n ! dimension of the vector v
      
! !INPUT/OUTPUT PARAMETERS:

      real(8), intent(inout) :: v(*)

      integer(4), optional, intent(inout) :: index(*)

! !LOCAL VARIABLES:

      integer(4) :: i,itmp,indtp,j
      real(8)    :: vtmp
      
! !SYSTEM ROUTINES:
      

! !REVISION HISTORY:
!
! Created 26th. Feb. 2004
!
!EOP
!BOC
      if(present(index))then
        do i=1,n
          index(i)=i
        enddo
      endif  
      do i=1,n-1
        vtmp=v(i)
        itmp=i
        do j=i+1,n
          if(v(j).lt.vtmp)then
            vtmp=v(j)
            itmp=j
          endif  
        enddo
         if(itmp.ne.i)then
          vtmp=v(i)
          v(i)=v(itmp)
          v(itmp)=vtmp
          if(present(index))then
            indtp=index(i)
            index(i)=index(itmp)
            index(itmp)=indtp
          endif  
        endif
      enddo
      end subroutine sortr
!EOC                
      end module order
