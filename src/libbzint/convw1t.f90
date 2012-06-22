
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: convw1t
!
! !INTERFACE: 
      subroutine convw1t(eo,eu,efer,weight)
      
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
      use tetra_internal, only: omgga, sgnfrq
 
! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in) :: eo(4) ! the band energies at k

      real(8), intent(in) :: eu(4) ! the band energies at k-q

      real(8), intent(in) :: efer   ! the fermi energy

! !OUTPUT PARAMETERS:

      real(8), intent(out) :: weight(4) ! the weight at each corner
       
! !LOCAL VARIABLES:

      integer(4) :: i,inod,ssg,inf
      integer(4) :: isub
      real(8), dimension(4,2) :: wtemp
      real(8), allocatable  :: nodtemp(:,:)
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
        if(sgnfrq.eq.4)then
          if(dabs(omgga).lt.1.0d-12) then
            call ksurf(e,ef,weight)
            do i=1,4
              weight(i)=weight(i)*(-3.1415926d0)
            enddo
          endif
        endif               

      case(4)     
        allocate(nodtemp(1:3,1:nnod))
        nodtemp(1:3,1:4)=intnodes(1:3,1:4)
        call generictetra(nodtemp,weight,1,inf)

      case(5) 
        allocate(nodtemp(1:3,1:nnod))
        do inod=1,5
          nodtemp(1:3,index(inod,1))=intnodes(1:3,inod)
        enddo
        call genericfunf(nodtemp,weight,1,inf)

      case(6)
        allocate(nodtemp(1:3,1:nnod))
        do inod=1,6
          nodtemp(1:3,index(inod,1))=intnodes(1:3,inod)
        enddo   
        call genericprism(nodtemp,weight,1,inf)

      case(7)
        allocate(nodtemp(1:3,1:5))
        do isub=1,2
          do inod=1,7
            if(index(inod,isub).ne.0)then
              nodtemp(1:3,index(inod,isub))=intnodes(1:3,inod)
            endif
          enddo ! inod
          call genericfunf(nodtemp,wtemp(1:4,isub),2,inf)
        enddo ! isub  
        do i=1,4
          weight(i)=wtemp(i,1)+wtemp(i,2)
        enddo
            
      case(8)
        allocate(nodtemp(1:3,1:6))
        do isub=1,2
          do inod=1,8
            if(index(inod,isub).ne.0) then
              nodtemp(1:3,index(inod,isub))=intnodes(1:3,inod)
            endif
          enddo ! inod
          call genericprism(nodtemp,wtemp(1:4,isub),4,inf)
        enddo ! isub  
        do i=1,4
          weight(i)=wtemp(i,1)+wtemp(i,2)
        enddo
       
      case default
         print*, nnod
         stop 'error in convw1t.f90(2)'
      end select

      if(allocated(nodtemp))deallocate(nodtemp)
      if(allocated(index))deallocate(index)
      
      end subroutine convw1t

!EOC        
