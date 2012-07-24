!BOP
!
! !ROUTINE: testmixfun
!
! !INTERFACE: 
      subroutine testmixfun
      
! !DESCRIPTION:
!
!This subroutine perform the test of the mixed functions, calculates the
!product of two eigenvectors directly and as a linear combination of mixed
!basis functions and writes both to disk for plotting.      
!
! !USES:

      use modmain
      use modgw

!
! !LOCAL VARIABLES:
!
      implicit none
      
      integer(4) :: ib1, ib2
      integer(4) :: ikp, jkp
      integer(4) :: ik, jk, iq
      
! !EXTERNAL ROUTINES: 

      external calcmpwipw
      external eptest
      external plotevecprod
      external plotevecmix
       
! !REVISION HISTORY:
!
! Created May 2006 by RGA
! Revisited 4.05.2011 by DIN
!
!EOP
!BOC
      call boxmsg(6,'-','TESTMIXFUN')
    
      ik=input%gw%iik
      jk=input%gw%jjk
            
      ikp=indkp(ik)
      jkp=indkp(jk)
      
      write(*,*) 'Parameters:'
      write(*,*) 'first k-point number (iik): ', input%gw%iik
      write(*,*) 'second k-point number (jjk): ', input%gw%jjk
      write(*,*) 'corresponding first irreducible k-point number (ikp): ', ikp
      write(*,*) 'corresponding second irreducible k-point number (jkp): ', jkp
      write(*,*) 'lower bound for band number (ibmin): ', input%gw%ibmin
      write(*,*) 'upper bound for band number (ibmax): ', input%gw%ibmax
      write(*,*) 'lower bound for band number (ibmin2): ', input%gw%ibmin2
      write(*,*) 'upper bound for band number (ibmax2): ', input%gw%ibmax2
      write(*,*) 'atom 1 (at1): ', input%gw%at1
      write(*,*) 'atom 2 (at2): ', input%gw%at2
      write(*,*)

!     find the corresponding q-vector id
      do iq=1,nqptnr
        if (kqid(ik,iq).eq.jk) exit
      enddo
      write(*,*) 'Corresponding Q-point: iq=', iq
      write(*,'(a,3f12.4)') 'vql=', vql(:,iq) 
!
      call diagsgi(iq)
!
      call calcmpwipw(iq)
!
      call eptest(ik,jk,iq)
!
      do ib1=input%gw%ibmin,input%gw%ibmax
        do ib2=input%gw%ibmin2,input%gw%ibmax2
!
          call plotevecprod(ik,jk,ib1,ib2,input%gw%at1,input%gw%at2)
!
          call plotevecmix(iq,ik,jk,ib1,ib2,input%gw%at1,input%gw%at2)
!
        enddo !ib2 
      enddo !ib1

      deallocate(sgi)
      deallocate(minmmat)
      deallocate(umix)
      deallocate(apwfr)
      deallocate(lofr)
      deallocate(locmixind)
      
      return
      end subroutine testmixfun

!EOC      
