
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)


!BOP
!
! !ROUTINE: convw_1k
!
! !INTERFACE:
subroutine convw_1k(iklib,efer,omeg,sigfreq,cweight)
!
! !DESCRIPTION:
!
!   This subroutine calculates the integration weight of each k-point for
!   all band pairs. Sigfreq distinguishes among the different kinds of weights.
!   sigfreq=1, normal q-dependent bulk integration.
!   sigfreq=2, weights for the Polarization with real frequencies
!   sigfreq=3, weights for the Polarization with imaginary frequencies.
!   sigfreq=4, it is for the q-dependent surface integration (the surface is
!   defined by e_jb-e_ib=omeg.
!   The convolution weight is calculated for all band combinations and one
!   {\bf k}-point.
!
! !USES:
  use order
  use tetra_internal
  implicit none
! !INPUT PARAMETERS:
  integer, intent(in) :: iklib
  real(8), intent(in)  :: efer
  real(8), intent(in)  :: omeg
  integer(4), intent(in)  :: sigfreq
! !OUTPUT PARAMETERS:
  real(8), intent(out) :: cweight(nband,nband) ! the weight
! !LOCAL VARIABLES:
  integer(4) :: itet,i,ib,jb,kin,kjn
  integer(4), dimension(4) :: ik1
  integer(4), dimension(4) :: ik2
  real(8) :: term,wwwt,tw
  real(8), dimension(4) :: ee1
  real(8), dimension(4) :: ee2
  real(8), dimension(4) :: w1t
  real(8), dimension(4) :: wcor
! !EXTERNAL ROUTINES:
  external intweight1t
  external convw1t

  external bloechlcor
! !SYSTEM ROUTINES:
  intrinsic maxval
  intrinsic minval
! !REVISION HISTORY:
!
!   Created: August 2008 by S. Sagmeister
!   Used {\tt convw.f90} as template for this routine.
!
!EOP
!BOC
  cweight=0.0d0
  omgga=omeg
  sgnfrq=sigfreq
  wwwt=0.0d0
  select case(sigfreq)
  case(1)
     ! normal q-dependent bulk integration
     do itet=1,ntet

        ! only specified k-point is considered
        if (all(tetcorn(:,itet).ne.iklib)) cycle

        tw=dble(tetweig(itet))
        do ib=1,nband
           do i=1,4
              ee1(i)=eband(ib,tetcorn(i,itet))
           end do
           if (maxval(ee1,dim=1).le.efer) then
              do jb=1,nband
                 do i=1,4
                    ee2(i)=eband(jb,tetcorn(i,tetln(itet)))
                 end do
                 if (minval(ee2,dim=1).gt.efer) then
                    do i=1,4
                       kin=tetcorn(i,itet)
                       kjn=tetcorn(i,tetln(itet))
                       ! if specified k-point matches
                       if (kin.eq.iklib) &
                            cweight(ib,jb)=cweight(ib,jb)+vt*tw/4.0d0
                    end do
                 else if (maxval(ee2,dim=1).gt.efer) then
                    do i=1,4
                       ee2(i)=efer-ee2(i)
                    end do
                    w1t(1:4)=0.0d0
                    wcor(1:4)=0.0d0
                    call sort(4,ee2,ik2)
                    call intweight1t(ee2,0.0d0,vt,w1t)
                    call bloechlcor(ee2,0.0d0,vt,wcor)
                    do i=1,4
                       term=w1t(i)+wcor(i)
                       kin=tetcorn(i,itet)                          !
                       kjn=tetcorn(ik2(i),tetln(itet))
                       if (kin.eq.iklib) &
                            cweight(ib,jb)=cweight(ib,jb)+term*tw
                    end do
                 end if
              end do
           else if (minval(ee1,dim=1).le.efer) then
              do jb=1,nband
                 do i=1,4
                    ee2(i)=eband(jb,tetcorn(i,tetln(itet)))
                 end do
                 if (minval(ee2,dim=1).gt.efer) then
                    w1t(1:4)=0.0d0
                    wcor(1:4)=0.0d0
                    call sort(4,ee1,ik1)
                    call intweight1t(ee1,efer,vt,w1t)
                    call bloechlcor(ee1,efer,vt,wcor)
                    do i=1,4
                       term=w1t(i)+wcor(i)
                       kin=tetcorn(ik1(i),itet)
                       kjn=tetcorn(i,tetln(itet))
                       if (kin.eq.iklib) &
                       cweight(ib,jb)=cweight(ib,jb)+term*tw
                    end do
                 else
                    w1t(1:4)=0.0d0
                    call convw1t(ee1,ee2,efer,w1t)
                    do i=1,4
                       kin=tetcorn(i,itet)
                       kjn=tetcorn(i,tetln(itet))
                       if (kin.eq.iklib) &
                            cweight(ib,jb)=cweight(ib,jb)+w1t(i)*vt*6*tw
                    end do
                 end if
              end do
           end if
        end do
     end do  ! itet
  case(2:4)
     ! for the q-dependent bulk integration for Polarization
     do itet=1,ntet

        ! only specified k-point is considered
        if (all(tetcorn(:,itet).ne.iklib)) cycle

        tw=dble(tetweig(itet))
        do ib=1,nband
           do i=1,4
              ee1(i)=eband(ib,tetcorn(i,itet))
           end do
           if (minval(ee1,dim=1).le.efer) then
              do jb=1,nband
                 do i=1,4
                    ee2(i)=eband(jb,tetcorn(i,tetln(itet)))
                 end do
                 if (maxval(ee2,dim=1).gt.efer) then
                    w1t(1:4)=0.0d0
                    call convw1t(ee1,ee2,efer,w1t)
                    do i=1,4
                       kin=tetcorn(i,itet)
                       kjn=tetcorn(i,tetln(itet))
                       if (kin.eq.iklib) &
                            cweight(ib,jb)=cweight(ib,jb)+w1t(i)*vt*6.0d0*tw
                    end do
                 end if
              end do
           end if
        end do
     end do
  case default
     continue
  end select
end subroutine convw_1k
!EOC
