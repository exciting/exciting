
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: convcorew
!
! !INTERFACE:
      subroutine convcorew(efer,omeg,sigfreq,cweight)
! 
! !DESCRIPTION:
!  
!   This subroutine calculates the integration weight of each k-point for
!   each band pair for the case when the occupied state is a core state.

!       
! !USES:
       
       use order
       
       use tetra_internal
       
       implicit none     
       
! !INPUT PARAMETERS:
 
       real(8), intent(in)  :: efer                       ! fermi energy
     
       real(8), intent(in)  :: omeg

       integer(4), intent(in)  :: sigfreq       

! !OUTPUT PARAMETERS:
       
       real(8), intent(out) :: cweight(nirkp,nband,ncore) ! the weight 
!                                                           of each 
!                                                           k-q point for
!                                                           each band and
!                                                           core state

!  
! !LOCAL VARIABLES:
 
       integer(4) :: itet,i,ib,jb,kin,ic
       integer(4), dimension(4) :: ik2
       real(8) :: term
       real(8), dimension(4) :: ee1
       real(8), dimension(4) :: ee2
       real(8), dimension(4) :: w1t
       real(8), dimension(4) :: wcor
       external  intweight1t
       
! !SYSTEM ROUTINES:
      intrinsic maxval
      intrinsic minval       
       
! !REVISION HISTORY:
!
!   Created:  17th. Jan 2005. by XZL 
!
!EOP
!BOC
 
      cweight(1:nirkp,1:nband,1:ncore)=0.0d0
      omgga=omeg
      sgnfrq=sigfreq

      select case (sigfreq)
      
      case(1)
        do itet=1,ntet
            do jb=1,nband
              do i=1,4
                ee2(i)=eband(jb,tetcorn(i,tetln(itet)))
              enddo

              if(minval(ee2,dim=1).gt.efer)then
                do i=1,4
                  kin=tetcorn(i,itet)
                  cweight(kin,jb,1:ncore)=cweight(kin,jb,1:ncore)+vt/4.0d0
                enddo
              else if(maxval(ee2,dim=1).gt.efer)then
                do i=1,4
                  ee2(i)=efer-ee2(i)
                enddo           
                w1t(1:4)=0.0d0
                wcor(1:4)=0.0d0
                call sort(4,ee2,ik2)
                call intweight1t(ee2,0.0d0,vt,w1t)
                call bloechlcor(ee2,0.0d0,vt,wcor)
                do i=1,4
                  term=w1t(i)+wcor(i)                           !
                  kin=tetcorn(ik2(i),itet)                             !
                  cweight(kin,jb,1:ncore)=cweight(kin,jb,1:ncore)+term
                enddo
              endif
            enddo
        enddo
          
      case(2:3)
        do itet=1,ntet
          do ic=1,ncore              
            do i=1,4
              ee1(i)=ecore(ic)
            enddo
            do jb=1,nband
               do i=1,4
                 ee2(i)=eband(2,tetcorn(i,tetln(itet)))
               enddo
               if(maxval(ee2,dim=1).gt.efer)then
!                 write(25,*)'rga: convw: calling bothpart1t with pars:'
!                 write(25,*)'rga: convw: ee1 = ',ee1
!                 write(25,*)'rga: convw: ee2 = ',ee2
!                 write(25,*)'rga: convw: efer = ',efer
                 w1t(1:4)=0.0d0
                 call bpartoc(ee1,ee2,efer,w1t)
!                 write(25,*)'rga: convw: bothpart1t returned:'
!                 write(25,*)'rga: convw: w1t = ',w1t
                 do i=1,4
                   kin=tetcorn(i,itet)
                   cweight(kin,jb,ic)=cweight(kin,ib,ic)+w1t(i)*vt*6.0d0
                 enddo
               else
                 continue
               endif
            enddo  
          enddo
        enddo
      end select
      
      return
      
      end subroutine convcorew

!EOC
