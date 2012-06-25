!BOP
!
! !ROUTINE: initgw
!
! !INTERFACE
subroutine testmbrotmat()

! !DESCRIPTION:
!
!
! 
! !USES:
      use modmain
      use modgw

! !LOCAL VARIABLES:
      implicit none
      integer(4) :: ik, iq, iqp
      integer(4) :: recl
      
! !EXTERNAL ROUTINES: 


! !REVISION HISTORY:
!
! Created, November 2011 by DIN
!
!EOP
!BOC
      iq=input%gw%iik
      iqp=indkp(iq) ! irreducible
      iq=idikp(iqp)
      
      write(*,*)
      write(*,*)'iq = ', iq, 'iqp = ', iqp
      write(*,*)

      write(*,*)
      call boxmsg(fgw,'-','q-dependent parameters')
      write(*,*)
!      
!     Calculate the interstitial mixed basis functions
!     
      call diagsgi(iq)
!
!     Set the size of the basis for the corresponding q-point
!      
      matsiz=locmatsiz+ngq(iq)
      write(*,101) iq, locmatsiz, ngq(iq), matsiz
!     
!     Calculate the transformation matrix between pw's and the mixed basis 
!
      call calcmpwipw(iq)
!      
!     Calculate the Minm matrix elements for the total (non-reduced) k-grid
!      
      recl=16*(locmatsiz*nstfv*ncg)
      open(39,file='mincmat.io',action='WRITE',form='UNFORMATTED', &
     &  access='DIRECT',recl=recl)

      recl=16*(locmatsiz*ncg*nstfv)
      open(40,file='micmmat.io',action='WRITE',form='UNFORMATTED', &
     &  access='DIRECT',recl=recl)
   
      recl=16*(matsiz*nstfv*nstsv)
      open(41,file='minmmat.io',action='WRITE',form='UNFORMATTED', &
     &  access='DIRECT',recl=recl)
   
      do ik = 1, nkptnr
        call expand_prods(ik,iq,0)
        write(39,rec=ik) mincmat
        write(40,rec=ik) micmmat
        write(41,rec=ik) minmmat
      end do
      close(39)
      close(40)
      close(41)
      
      call checkmbrot(iq)

  101 format(10x,'Data for q-point nr.:',i4,//,10x,'Mixed basis:',/,10x, &
     &           'Number of atomic basis functions:       ',i4,/,10x,   &
     &           'Number of interstitial basis functions: ',i4,/,10x,   &
     &           'Total number of basis functions:        ',i4,/)
      
      return
end subroutine
!EOC
