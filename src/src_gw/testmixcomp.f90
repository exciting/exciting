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
      integer(4) :: ikp, jkp
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
      call boxmsg(6,'-','TESTMIXCOMP')

      ik=input%gw%iik
      jk=input%gw%jjk
      
      ikp=indkp(ik)
      jkp=indkp(jk)
      
      write(*,*) 'Parameters:'
      write(*,*) 'first k-point number (iik): ', ik
      write(*,*) 'second k-point number (jjk): ', jk
      write(*,*) 'corresponding first irreducible k-point number (ikp): ', ikp
      write(*,*) 'corresponding second irreducible k-point number (jkp): ', jkp
      write(*,*) 'lower bound for band number (ibmin): ', input%gw%ibmin
      write(*,*) 'upper bound for band number (ibmax): ', input%gw%ibmax
      write(*,*) 'lower bound for band number (ibmin2): ', input%gw%ibmin2
      write(*,*) 'upper bound for band number (ibmax2): ', input%gw%ibmax2
      write(*,*)

!------------------------------------------------------------------------
      
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

      call cpu_time(tt(1))
      call diagsgi(iq)
      call cpu_time(tt(2))
      call write_cputime(fgw,tt(2)-tt(1),'DIAGSGI')

      call cpu_time(tt(1))
      call calcmpwipw(iq)
      call cpu_time(tt(2))
      call write_cputime(fgw,tt(2)-tt(1),'CALCMPIPW')

      call cpu_time(tt(1))
      call eptest(ik,jk,iq)
      call cpu_time(tt(2))
      call write_cputime(fgw,tt(2)-tt(1),'EPTEST') 
      
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
      call write_cputime(fgw,tt(2)-tt(1),'INTEVECPP')

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
      call write_cputime(fgw,tt(2)-tt(1),'INTEVECPM')

      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)

      return
      end subroutine testmixcomp
!EOC
