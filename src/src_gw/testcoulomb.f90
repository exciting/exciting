!BOP
!
! !ROUTINE: testcoulomb
!
! !INTERFACE:
      subroutine testcoulomb

! !DESCRIPTION:
!
! This subroutine performs the test of the bare coulomb potential,
! calculates it with the mixed basis and planewaves and compares both.
!
! !USES:

     use modinput
     use modmain
     use modgw
           
! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: iq

! !EXTERNAL ROUTINES: 

      real(8) calceta
      real(8) gcutoff
      real(8) rcutoff

! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
! Revisited 23.05.2011 by DIN
!
!EOP
!BOC
      open(unit=93,file='barcmix.evec',form='formatted',status='unknown')
      open(unit=94,file='barcpw.evec',form='formatted',status='unknown')
      open(unit=95,file='barcr.mat',form='formatted',status='unknown')
      open(unit=96,file='barci.mat',form='formatted',status='unknown')
      open(unit=97,file='barc.eval',form='formatted',status='unknown')
      open(unit=99,file='TESTBARC.OUT',form='formatted',status='unknown')

      do iq = input%gw%iik, input%gw%jjk

        write(93,*)'#----------------------------------------------'
        write(93,*)'#              iq = ', iq
        write(93,*)'#----------------------------------------------'
        write(94,*)'#----------------------------------------------'
        write(94,*)'#              iq = ', iq
        write(94,*)'#----------------------------------------------'
        write(95,*)'#----------------------------------------------'
        write(95,*)'#              iq = ', iq
        write(95,*)'#----------------------------------------------'
        write(96,*)'#----------------------------------------------'
        write(96,*)'#              iq = ', iq
        write(96,*)'#----------------------------------------------'
        write(97,*)'#----------------------------------------------'
        write(97,*)'#              iq = ', iq
        write(97,*)'#----------------------------------------------'
        
        Gamma=gammapoint(iq)
        
        call diagsgi(iq)
        
        matsiz=locmatsiz+ngq(iq)
        write(*,101)iq,locmatsiz,ngq(iq),matsiz
        
        call calcmpwipw(iq)
        
        if (Gamma) then
          call calcwmix0
        end if
        
        call testbarc(iq)
        
        deallocate(mpwmix)
        deallocate(mpwipw)
        deallocate(barc)
        deallocate(sqbarc)

      enddo  
      close(93)
      close(94)
      close(95)
      close(96)
      close(97)
      close(99)
      close(52)

  101 format(10x,'Data for q-point nr.:',i4,//,10x,'Mixed basis:',/,10x, &
     &           'Number of atomic basis functions:       ',i4,/,10x,   &
     &           'Number of interstitial basis functions: ',i4,/,10x,   &
     &           'Total number of basis functions:        ',i4,/)
      
      return
      end subroutine testcoulomb

!EOC      
