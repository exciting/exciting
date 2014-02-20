!BOP
!
! !ROUTINE: calcmpwipw
!
! !INTERFACE:
      subroutine calcmpwipw(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the matrix elements between PW's and
! orthogonalized IPW's.
!
! !USES:

      use modmain
      use modgw
      use modmpi, only: rank
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: iq

! !LOCAL VARIABLES:

      integer(4) :: ipw   ! (Counter): runs over plane waves
      integer(4) :: jpw   ! (Coutner): runs over IPW's
      integer(4) :: ig

      integer(4), dimension(3) :: igv ! integer coodintates of G-G'

      complex(8), allocatable  :: tmat(:,:)
      
      real(8) :: tstart, tend
 
! !EXTERNAL ROUTINES: 

      external zgemm
!
! !REVISION HISTORY:
! 
! Created: May 2006 by RGA
! Revisited 10.05.2011 by DIN

!EOP
!BOC

      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

      if(allocated(mpwipw))deallocate(mpwipw)
      allocate(mpwipw(ngq(iq),ngbarc(iq)))
      allocate(tmat(ngq(iq),ngbarc(iq)))
      mpwipw=zzero
!
!     Calculate the integral between pw's and IPW's
!
      tmat(:,:)=zzero
      do ipw=1,ngbarc(iq)
        do jpw=1,ngq(iq)
          igv(:)=ivg(:,igqigb(ipw,iq))-ivg(:,igqig(jpw,iq))
          if((igv(1).ge.intgv(1,1)).and.(igv(1).le.intgv(1,2)).and.       &
         &   (igv(2).ge.intgv(2,1)).and.(igv(2).le.intgv(2,2)).and.       &
         &   (igv(3).ge.intgv(3,1)).and.(igv(3).le.intgv(3,2)))        then
            ig=ivgig(igv(1),igv(2),igv(3))
            !tmat(jpw,ipw)=conjg(cfunig(ig))
            tmat(jpw,ipw)=ipwint(ig)
          else
            tmat(jpw,ipw)=zzero
          end if
        enddo ! jpw  
      enddo ! ipw
      
!     Calculate the overlap PW-IPW integral
      call zgemm('t','n',ngq(iq),ngbarc(iq),ngq(iq),zone,sgi,ngq(iq), &
     &           tmat,ngq(iq),zzero,mpwipw,ngq(iq))
        
      if (debug) then
           write(55,*) 'CALCMPWIPW, iq = ', iq
           write(55,*) 'CALCMPWIPW, ngq = ', ngq(iq)
           write(55,*) 'CALCMPWIPW, ngbarc = ', ngbarc(iq)
           write(55,*) "### mpwipw for iq=",iq 
           write(55,*)   
           write(55,*) 'jpw  ipw  ivg(ipw)  ivg(jpw)  cfunig &
          &             mpwipw(jpw,ipw)'

           do ipw=1,ngbarc(iq),ngbarc(iq)/10
             do jpw=1,ngq(iq),ngq(iq)/10
               igv(:)=ivg(:,igqigb(ipw,iq))-ivg(:,igqig(jpw,iq))
               write(55,'(8i5,4f12.6)') jpw, ipw,               &
          &         ivg(:,igqigb(ipw,iq)),ivg(:,igqig(jpw,iq)), &
          &         tmat(jpw,ipw),mpwipw(jpw,ipw)
             enddo 
           enddo 
           
           write(55,*) "### sgi ###"
           do jpw=1,ngq(iq),ngq(iq)/10
             do ipw=1,ngq(iq),ngq(iq)/10
               write(55,'(2i5,4f12.6)') ipw,jpw,sgi(ipw,jpw)
             enddo 
           enddo  
      end if ! debug
      
      deallocate(tmat)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if (rank==0) call write_cputime(fgw,tend-tstart, 'CALCMPWIPW')

      return
      end subroutine calcmpwipw
!EOC
