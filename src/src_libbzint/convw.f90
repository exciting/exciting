!BOP
!
! !ROUTINE: convw
!
! !INTERFACE:
      subroutine convw(efer,omeg,sigfreq,cwt)
! 
! !DESCRIPTION:
!  
!   This subroutine calculates the integration weight of each k-point for
!   all band pairs. Sigfreq distinguishes among the different kinds of weights.
!   sigfreq=1, normal q-dependent bulk integration. 
!   sigfreq=2, weights for the Polarization with real frequencies
!   sigfreq=3, weights for the Polarization with imaginary frequencies.
!   sigfreq=4, it is for the q-dependent surface integratio (the surface is defined by
!   e_jb-e_ib=omeg.

!       
! !USES:
       
       use order
       
       use tetra_internal
       
       implicit none     
       
! !INPUT PARAMETERS:
       real(8), intent(in)  :: efer                       ! fermi energy
       real(8), intent(in)  :: omeg
       integer, intent(in)  :: sigfreq       
! !OUTPUT PARAMETERS:
       real(8), intent(out) :: cwt(nband,nband,nirkp) ! the weight of each k-point for each band 
 
!  
! !LOCAL VARIABLES:
 
      integer :: it,i,ib,jb,kin
      integer :: ik1(4),ik2(4) 
      real(8) :: term, wwwt, tw
      real(8) :: ee1(4),ee2(4),w1t(4),wcor(4) 

! !EXTERNAL ROUTINES:
      external intweight1t
      external convw1t
      external convw1tsurf
      external bloechlcor

! !SYSTEM ROUTINES:
      intrinsic maxval
      intrinsic minval       
       
! !REVISION HISTORY:
!
!   Created:  3th. March 2004. by RGA 
!
!EOP
!BOC
 
      cwt=0.0d0
      omgga=omeg
      sgnfrq=sigfreq
      wwwt=0.0d0
      select case (sigfreq)
      
      case(1)      ! normal q-dependent bulk integration
      do it=1,ntet
        tw=dble(tetweig(it))
        do ib=1,nband
          do i=1,4
            ee1(i)=eband(ib,tetcorn(i,it))
          enddo
          if(maxval(ee1).le.efer)then
            do jb=1,nband          
              do i=1,4
                ee2(i)=eband(jb,tetcorn(i,tetln(it)))
              enddo

              if(minval(ee2).gt.efer)then
                do i=1,4
                  kin=tetcorn(i,it)
                  cwt(ib,jb,kin)=cwt(ib,jb,kin)+vt*tw/4.0d0
                enddo
              else if(maxval(ee2).gt.efer)then
                do i=1,4
                  ee2(i)=efer-ee2(i)
                enddo           
                w1t(1:4)=0.0d0
                wcor(1:4)=0.0d0
                call sort(4,ee2,ik2)
                call intweight1t(ee2,0.0d0,vt,w1t)
                call bloechlcor(ee2,0.0d0,vt,wcor)
                do i=1,4
                  term=w1t(i)+wcor(i)
                  kin=tetcorn(i,it)                          !
                  cwt(ib,jb,kin)=cwt(ib,jb,kin)+term*tw
                enddo
              endif
            enddo
          
          else if(minval(ee1).le.efer)then
            do jb=1,nband
              do i=1,4
                ee2(i)=eband(jb,tetcorn(i,tetln(it)))
              enddo
              if(minval(ee2).gt.efer)then
                w1t(1:4)=0.0d0
                wcor(1:4)=0.0d0
                call sort(4,ee1,ik1)
                call intweight1t(ee1,efer,vt,w1t)
                call bloechlcor(ee1,efer,vt,wcor)
                do i=1,4
                  term=w1t(i)+wcor(i)
                  kin=tetcorn(ik1(i),it)
                  cwt(ib,jb,kin)=cwt(ib,jb,kin)+term*tw
                enddo
              else 
                  w1t(1:4)=0.0d0
                  call convw1t(ee1,ee2,efer,w1t)
                  do i=1,4
                    kin=tetcorn(i,it)
                    cwt(ib,jb,kin)=cwt(ib,jb,kin)+w1t(i)*vt*6*tw
                  enddo
              
              endif
            enddo
          endif
        enddo
      enddo  ! it

      case(2:4)     ! for the q-dependent bulk integration for Polarization                      
        do it=1,ntet
          tw=dble(tetweig(it))
          do ib=1,nband  ! occupied band 
            do i=1,4
              ee1(i)=eband(ib,tetcorn(i,it))
            enddo
            if(minval(ee1).gt.efer) cycle 

            do jb=1,nband  ! unoccupied band 
              do i=1,4
                ee2(i)=eband(jb,tetcorn(i,tetln(it)))
              enddo
              if(maxval(ee2).le.efer) cycle 

              w1t(1:4)=0.0d0
              call convw1t(ee1,ee2,efer,w1t)
              do i=1,4
                kin=tetcorn(i,it)
                cwt(ib,jb,kin)=cwt(ib,jb,kin)+w1t(i)*vt*6.0d0*tw
              enddo
            enddo
          enddo
        enddo
      
      case default
        continue
      end select
      
      return
      
      end subroutine convw

!EOC
