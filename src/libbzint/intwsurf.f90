
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: intwsurf
!
! !INTERFACE:
       subroutine intwsurf(omeg,iweight)
 
! !DESCRIPTION:
!
!   This subroutine calculates the surface integration weight of each k-point for
!   each band

! !USES:
       
       use order
       use tetra_internal
       
       implicit none     

! !INPUT PARAMETERS:
 
       real(8), intent(in)  :: omeg                 ! fermi energy
       
! !OUTPUT PARAMETERS:

       real(8), intent(out) :: iweight(nband,nirkp) ! the weight of each
!                                                     k-point for each 
!                                                     band

!  
! !LOCAL VARIABLES:
 
       integer(4) :: itet,i,ib,kin
       integer(4), dimension(4) :: ik
       real(8) :: term
       real(8), dimension(4) :: ee
       real(8), dimension(4) :: w1t

! !EXTERNAL ROUTINES:

       external  intweight1t

! !REVISION HISTORY:
!
!   Created:  10th. Jan 2005. by XZL 
! 
!EOP
!BOC
      do i=1,nirkp
        do ib=1,nband
          iweight(ib,i)=0.0d0
        enddo
      enddo
   
      do itet=1,ntet
        do ib=1,nband 
          do i=1,4
            ee(i)=eband(ib,tetcorn(i,itet))
          enddo
          call sort(4,ee,ik)
          w1t(1:4)=0.0d0
          if((ee(1).le.omeg).and.(omeg.le.ee(4))) then
            call ksurf(ee,omeg,w1t)     
            do i=1,4
             term=w1t(i)*tetweig(itet)
             kin=tetcorn(ik(i),itet)
             iweight(ib,kin)=iweight(ib,kin)+term*6.0d0*vt
            enddo
          else
            continue
          endif
        enddo
      enddo

      return
      
      end subroutine intwsurf
      
!EOC
