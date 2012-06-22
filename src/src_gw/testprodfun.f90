!BOP
!
! !ROUTINE: testprodfun
!
! !INTERFACE: 
      subroutine testprodfun
      
! !DESCRIPTION:
!
!This subroutine perform the test of the mixed functions, calculates the
!product of two eigenvectors directly and as a linear combination of mixed
!basis functions and writes both to disk for plotting.      
!
! !USES:
      use modinput
!
! !LOCAL VARIABLES:
!
      implicit none
      integer(4) :: ik, jk
      integer(4) :: ib1, ib2

! !EXTERNAL ROUTINES: 

      external plotevecprod
       
! !REVISION HISTORY:
!
! Created 16.09.05 by RGA
! Revisited 03.05.2011 by DIN
!
!EOP
!BOC 

      call boxmsg(6,'-','TESTPRODFUN')
      
      write(*,*) 'Parameters:'
      write(*,*) 'first k-point number (iik): ', input%gw%iik
      write(*,*) 'second k-point number (jjkk): ', input%gw%jjk
      write(*,*) 'lower bound for band number (ibmin): ', input%gw%ibmin
      write(*,*) 'upper bound for band number (ibmax): ', input%gw%ibmax
      write(*,*) 'lower bound for band number (ibmin2): ', input%gw%ibmin2
      write(*,*) 'upper bound for band number (ibmax2): ', input%gw%ibmax2
      write(*,*) 'atom 1 (at1): ', input%gw%at1
      write(*,*) 'atom 2 (at2): ', input%gw%at2
      write(*,*)
 
      ik=input%gw%iik
      jk=input%gw%jjk
      do ib1=input%gw%ibmin,input%gw%ibmax
        do ib2=input%gw%ibmin2,input%gw%ibmax2
          call plotevecprod(ik,jk,ib1,ib2,input%gw%at1,input%gw%at2)
        enddo !ib2 
      enddo !ib1
      
      return
      end subroutine testprodfun

!EOC      
