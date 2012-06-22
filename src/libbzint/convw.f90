
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: convw
!
! !INTERFACE:
      subroutine convw(efer,omeg,sigfreq,cweight)
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

       integer(4), intent(in)  :: sigfreq

! !OUTPUT PARAMETERS:

       real(8), intent(out) :: cweight(nband,nband,nirkp) ! the weight
!                                                           of each
!                                                           k-point for
!                                                           each band

!
! !LOCAL VARIABLES:

      integer(4) :: itet,i,ib,jb,kin,kjn
      integer(4), dimension(4) :: ik1
      integer(4), dimension(4) :: ik2
      real(8) :: term, wwwt, tw
      real(8), dimension(4) :: ee1
      real(8), dimension(4) :: ee2
      real(8), dimension(4) :: w1t
      real(8), dimension(4) :: wcor

! !EXTERNAL ROUTINES:

      external intweight1t
      external convw1t

!      external sort
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

      cweight=0.0d0
      omgga=omeg
      sgnfrq=sigfreq
      wwwt=0.0d0
      select case (sigfreq)

      case(1)      ! normal q-dependent bulk integration
      do itet=1,ntet
        tw=dble(tetweig(itet))
        do ib=1,nband
          do i=1,4
            ee1(i)=eband(ib,tetcorn(i,itet))
          enddo
          if(maxval(ee1,dim=1).le.efer)then
            do jb=1,nband
              do i=1,4
                ee2(i)=eband(jb,tetcorn(i,tetln(itet)))
              enddo

              if(minval(ee2,dim=1).gt.efer)then
                do i=1,4
                  kin=tetcorn(i,itet)
                  kjn=tetcorn(i,tetln(itet))
                  cweight(ib,jb,kin)=cweight(ib,jb,kin)+vt*tw/4.0d0
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
                  term=w1t(i)+wcor(i)
                  kin=tetcorn(i,itet)                          !
                  kjn=tetcorn(ik2(i),tetln(itet))
                  cweight(ib,jb,kin)=cweight(ib,jb,kin)+term*tw
                enddo
              endif
            enddo

          else if(minval(ee1,dim=1).le.efer)then
            do jb=1,nband
              do i=1,4
                ee2(i)=eband(jb,tetcorn(i,tetln(itet)))
              enddo
              if(minval(ee2,dim=1).gt.efer)then
                w1t(1:4)=0.0d0
                wcor(1:4)=0.0d0
                call sort(4,ee1,ik1)
                call intweight1t(ee1,efer,vt,w1t)
                call bloechlcor(ee1,efer,vt,wcor)
                do i=1,4
                  term=w1t(i)+wcor(i)
                  kin=tetcorn(ik1(i),itet)
                  kjn=tetcorn(i,tetln(itet))
                  cweight(ib,jb,kin)=cweight(ib,jb,kin)+term*tw
                enddo
              else
!                  write(25,*)'rga: convw: calling bothpart1t with pars:'
!                  write(25,*)'rga: convw: ee1 = ',ee1
!                  write(25,*)'rga: convw: ee2 = ',ee2
!                  write(25,*)'rga: convw: efer = ',efer
                  w1t(1:4)=0.0d0
                  call convw1t(ee1,ee2,efer,w1t)
!                  write(25,*)'rga: convw: bothpart1t returned:'
!                  write(25,*)'rga: convw: w1t = ',w1t
                  do i=1,4
                    kin=tetcorn(i,itet)
                    kjn=tetcorn(i,tetln(itet))
                    cweight(ib,jb,kin)=cweight(ib,jb,kin)+w1t(i)*vt*6*tw
                  enddo

              endif
            enddo
          endif
        enddo
      enddo  ! itet


      case(2:4)     ! for the q-dependent bulk integration for Polarization
      do itet=1,ntet
        tw=dble(tetweig(itet))
        do ib=1,nband
          do i=1,4
            ee1(i)=eband(ib,tetcorn(i,itet))
          enddo
          if(minval(ee1,dim=1).le.efer)then
            do jb=1,nband
              do i=1,4
                ee2(i)=eband(jb,tetcorn(i,tetln(itet)))
              enddo
              if(maxval(ee2,dim=1).gt.efer)then
                w1t(1:4)=0.0d0
                call convw1t(ee1,ee2,efer,w1t)
                do i=1,4
                 kin=tetcorn(i,itet)
                 kjn=tetcorn(i,tetln(itet))
                 cweight(ib,jb,kin)=cweight(ib,jb,kin)+w1t(i)*vt*6.0d0*tw
                enddo
              endif
            enddo
          endif
        enddo
      enddo

      case default
        continue
      end select

      return

      end subroutine convw

!EOC
