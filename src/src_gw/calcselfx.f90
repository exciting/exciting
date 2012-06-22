!BOP
!
! !ROUTINE: addtovx
!
! !INTERFACE: 
      subroutine calcselfx(ikp,iqp)
      
! !DESCRIPTION:
!
!This subroutine calculates the (k,q) exchange term of the 
!self energy
!
! !USES:
     
      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: ikp
      integer(4), intent(in) :: iqp
      
! !LOCAL VARIABLES:            

      integer(4) :: icg     ! (Counter) Runs over core states.
      integer(4) :: ie1,ie2 ! (Counter) Runs over bands.
      integer(4) :: m,m1,m2

      real(8) :: tstart, tend, tcore
      real(8) :: sxs2
      
      complex(8) :: mvm     ! Sum_ij{M^i*W_ij(w)*conjg(M^j)}
      complex(8), allocatable :: sxqval(:),sxqcor(:)      
      complex(8) :: sum

! !INTRINSIC ROUTINES: 

      intrinsic abs
      intrinsic sign
      intrinsic cmplx
      intrinsic conjg
      intrinsic cpu_time

! !EXTERNAL ROUTINES: 

      complex(8), external :: zdotc
      external zgemv
 
! !REVISION HISTORY:
!
! Created 23.06.05 by RGA.
! Revisited June 2011 by DIN
!
!EOP
!BOC      
      call cpu_time(tstart)
     
      sxs2=-4.d0*pi*vi

!-------------------------------- 
!       Valence contribution
!-------------------------------- 
      
      allocate(sxqval(ibgw:nbgw))
      sxqval(ibgw:nbgw)=zzero
!
!     Summation over occupied bands
!
      do ie1 = ibgw, nbgw
         m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
         m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
         sum=zzero
         do m = m1, m2
           do ie2 = 1, maxoccband
             mvm=zdotc(mbsiz,minmmat(1:mbsiz,m,ie2),1,minmmat(1:mbsiz,m,ie2),1)
             sum=sum+mvm
           enddo ! ie2
         end do ! m           
         sxqval(ie1)=sxqval(ie1)-sum*wkpq(iqp,ikp)/dble(m2-m1+1)
         if ((iqp.eq.1).and.(ie1.le.maxoccband)) then
            sxqval(ie1)=sxqval(ie1)+sxs2*singc2
         end if
      enddo ! ie1

!-------------------------------- 
!       Core electron contribution
!-------------------------------- 

      call cpu_time(tcore)

      if (wcore) then

        allocate(sxqcor(ibgw:nbgw))
        sxqcor(ibgw:nbgw)=zzero
!
!       Summation over occupied bands
! 
        do ie1 = ibgw, nbgw
           sum=zzero
           m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
           m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
           do m = m1, m2
             do icg = 1, ncg
               mvm=zdotc(mbsiz,mincmat(1:mbsiz,m,icg),1,mincmat(1:mbsiz,m,icg),1)
               sum=sum+mvm
             enddo ! icg
           end do ! m
           sxqcor(ie1)=sxqcor(ie1)-sum*wkpq(iqp,ikp)/dble(m2-m1+1)
        enddo ! ie1
 
      endif ! wcore
      
!----------------------------------------
!     Sum up the contributions
!----------------------------------------

      write(96,*)'-------- CALCSELFX -------------, ikp = ', ikp, '    iqp = ', iqp
      if (wcore) then
        write(96,*)'# band nr.       core           valence        selfex'
        do ie1 = ibgw, nbgw
          selfex(ie1,ikp)=selfex(ie1,ikp)+sxqval(ie1)+sxqcor(ie1)
          write(96,10)ie1,sxqcor(ie1),sxqval(ie1),selfex(ie1,ikp)
        enddo ! ie1
        write(96,*)
        deallocate(sxqval)
        deallocate(sxqcor)
      else
        write(96,*)'# band nr.       valence        selfex'  
        do ie1 = ibgw, nbgw
           selfex(ie1,ikp)=selfex(ie1,ikp)+sxqval(ie1)
           if((iqp.eq.1).and.(ie1.le.maxoccband))then
             selfex(ie1,ikp)=selfex(ie1,ikp)+sxs2*singc2
           end if
           write(96,11)ie1,sxqval(ie1),selfex(ie1,ikp)
        enddo ! ie1
        write(96,*)
        deallocate(sxqval)
      endif

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
      call write_cputime(fgw,tend-tstart, 'CALCSELFX')
      call write_cputime(fgw,tcore-tstart,'CALCSELFX: core')
      call write_cputime(fgw,tend-tcore,  'CALCSELFX: valence')

   10 format(i4,2f18.10,'  | ',2f18.10,' || ',2f18.10)
   11 format(i4,2f18.10,' || ',2f18.10)
      
      return 
      end subroutine calcselfx
!EOC        
              
                 
                    
                    
                               
                        
