
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: intw 
!
! !INTERFACE:
       subroutine intw(efer,iweight)
!
! !DESCRIPTION:
!
!   This subroutine calculates the integration weight of each k-point for
!   each band

! !USES:
       
       use order
       use tetra_internal
       
       implicit none     

! !INPUT PARAMETERS:
 
       real(8), intent(in)  :: efer                 ! fermi energy
       
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
       real(8), dimension(4) :: wcor
       external  intweight1t

! !REVISION HISTORY:
!
!   Created:  3th. March 2004. by RGA 
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
          call intweight1t(ee,efer,vt,w1t)
          call bloechlcor(ee,efer,vt,wcor)
          
          do i=1,4
            term=(w1t(i)+wcor(i))*tetweig(itet)
            kin= tetcorn(ik(i),itet)
            iweight(ib,kin)=iweight(ib,kin)+term
          enddo
        enddo
      enddo
      
      return
      
      end subroutine intw
      
!EOC
