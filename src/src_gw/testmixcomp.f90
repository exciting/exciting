!BOP
!
! !ROUTINE: testmixcomp
!
! !INTERFACE:
      subroutine testmixcomp
      
! !DESCRIPTION:
!
! This subroutine tests the completeness of the mixed basis
!      
! !USES:
      
      use modinput
      use modmain
      use modgw
      
! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: ik, jk, iq     
      integer(4) :: ib1, ib2
      integer(4) :: ml,ngr

      real(8) ::  tt(2)
            
      character(len=64) :: fnat
      character(len=64) :: fnint
      character(len=64) :: fntot

! !INTRINSIC ROUTINES:
      
      intrinsic cpu_time     

! !REVISION HISTORY:
!     
! Created 16. Sept 2005 by RGA
!
!EOP
!BOC

!------------------------------------------------------------------------

      ik=input%gw%iik
      jk=input%gw%jjk
      
  71  format('intevppat-',i4,'-',i4,'.out')
      write(fnat,71) ik, jk
      call str_strip(fnat)

  72  format('intevppi-',i4,'-',i4,'.out')
      write(fnint,72) ik, jk
      call str_strip(fnint)
   
  73  format('intevpptot-',i4,'-',i4,'.out')
      write(fntot,73) ik, jk
      call str_strip(fntot)

      open(unit=71,file=fnat,status='unknown')
      open(unit=72,file=fnint,status='unknown')
      open(unit=73,file=fntot,status='unknown')

!------------------------------------------------------------------------

  74  format('intevpmat-',i4,'-',i4,'.out')
      write(fnat,74) ik, jk
      call str_strip(fnat)

  75  format('intevpmi-',i4,'-',i4,'.out')
      write(fnint,75) ik, jk
      call str_strip(fnint)
   
  76  format('intevpmtot-',i4,'-',i4,'.out')
      write(fntot,76) ik, jk
      call str_strip(fntot)

      open(unit=74,file=fnat,status='unknown')
      write(74,*)'ib1  ib2  ias    epint   emint  abserr    relerr'
      open(unit=75,file=fnint,status='unknown')
      write(75,*)'ib1  ib2  epii   emii           abserr    relerr'
      open(unit=76,file=fntot,status='unknown')
      write(76,*)'ib1  ib2  intep  intmix         abserr    relerr'

!------------------------------------------------------------------------

!     find the corresponding q-vector id
      do iq=1,nqptnr
        if (kqid(ik,iq).eq.jk) exit
      enddo  

      open(98,file='TESTMIXCOMP.OUT',action='write')

      call cpu_time(tt(1))
      call diagsgi(iq)
      call cpu_time(tt(2))
      write(98,101)tt(2)-tt(1)

      call cpu_time(tt(1))
      call calcmpwipw(iq)
      call cpu_time(tt(2))
      write(98,102)tt(2)-tt(1)

      call cpu_time(tt(1))
      call eptest(ik,jk,iq)
      call cpu_time(tt(2))
      write(98,103)tt(2)-tt(1)
      
! ----------------------------------------------------------
      ml=input%groundstate%lmaxapw
      ngr=16*ml*ml/3
      call prep_ang_int(ml,ngr)
      
      call cpu_time(tt(1))
      do ib1=input%gw%ibmin,input%gw%ibmax
        do ib2=input%gw%ibmin2,input%gw%ibmax2
          call intevecpp(ik,jk,ib1,ib2)
        enddo  
      enddo  
      call cpu_time(tt(2))
      write(98,104)tt(2)-tt(1)

      rewind(71)
      rewind(72)
      rewind(73)

! ----------------------------------------------------------
      call cpu_time(tt(1))
      do ib1=input%gw%ibmin,input%gw%ibmax
        do ib2=input%gw%ibmin2,input%gw%ibmax2
          call intevecpm(ib1,ib2)
        enddo  
      enddo  
      call cpu_time(tt(2))
      write(98,105)tt(2)-tt(1)

      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)
  101 format(' diagsgi    -->',f14.5,' seconds')
  102 format(' calcmpwipw -->',f14.5,' seconds')
  103 format(' eptest     -->',f14.5,' seconds')
  104 format(' integral prods -->',f14.5,' seconds')
  105 format(' integral mix   -->',f14.5,' seconds')
      return
      
      end subroutine testmixcomp
            
!EOC
