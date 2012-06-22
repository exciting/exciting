!BOP
!
! !ROUTINE: addtovx
!
! !INTERFACE: 
      subroutine selfxsing(ikp)
      
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
      
! !LOCAL VARIABLES:            

      integer(4) :: icg     ! (Counter) Runs over core states.
      integer(4) :: ie1,ie2,ie12 ! (Counter) Runs over bands.
      integer(4) :: m,m1,m2
      integer(4) :: dimtk
     
      real(8) ::    tstart, tend, tcore
      
      complex(8) :: mvm     ! Sum_ij{M^i*W_ij(w)*conjg(M^j)}
      complex(8), allocatable :: vect1(:),vect2(:)
      complex(8), allocatable :: mmat1(:,:),mmat2(:,:)
      complex(8), allocatable :: lbarc(:,:)
      
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
     
!-------------------------------- 
!       Valence contribution
!-------------------------------- 

      dimtk=nbandsgw*nstfv
      allocate(mmat1(1:matsiz,1:dimtk))
      allocate(mmat2(1:matsiz,1:dimtk))
        
      ie12=0
      do ie1 = ibgw, nbgw
        do ie2 = 1, maxoccband
          ie12=ie12+1
          mmat1(1:matsiz,ie12)=minmmat(1:matsiz,ie1,ie2)
        end do
      end do  
!
!       v_{ij}*M^i_{nm}
!
      allocate(lbarc(1:matsiz,1:matsiz))

      lbarc(1:matsiz,1:matsiz)=barcs(1:matsiz,1:matsiz)
      call zhemm('l','u',matsiz,dimtk,-zone,lbarc,matsiz,      &
     &  mmat1,matsiz,zzero,mmat2,matsiz)

      deallocate(lbarc)
!
!     Summation over occupied bands
!
      allocate(vect1(1:matsiz))
      allocate(vect2(1:matsiz))
      do ie1 = ibgw, nbgw
         sum=zzero      
         m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
         m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
         do m = m1, m2
           do ie2 = 1, maxoccband
             ie12=(m-ibgw)*maxoccband+ie2
             vect1(1:matsiz)=mmat1(1:matsiz,ie12)
             vect2(1:matsiz)=mmat2(1:matsiz,ie12)
             mvm=zdotc(matsiz,vect1,1,vect2,1)
             sum=sum+mvm
           enddo ! ie2
         end do ! m           
         selfxs2(ie1,ikp)=selfxs2(ie1,ikp)+sum/(m2-m1+1)
      enddo ! ie1
      deallocate(vect1,vect2)
      deallocate(mmat1,mmat2)

!-------------------------------- 
!       Core electron contribution
!-------------------------------- 

      call cpu_time(tcore)

      if (wcore) then

        dimtk=nbandsgw*ncg
        allocate(mmat1(1:locmatsiz,1:dimtk))
        allocate(mmat2(1:locmatsiz,1:dimtk))       
       
        ie12=0
        do ie1 = ibgw, nbgw
           do icg = 1, ncg
             ie12=ie12+1
             mmat1(1:locmatsiz,ie12)=mincmat(1:locmatsiz,ie1,icg)
           enddo ! icg
        enddo ! ie1
!
!       v_{ij}*M^i_{nm}
!
        allocate(lbarc(1:locmatsiz,1:locmatsiz))

        lbarc(1:locmatsiz,1:locmatsiz)=barcs(1:locmatsiz,1:locmatsiz)
        call zhemm('l','u',locmatsiz,dimtk,-zone,lbarc,locmatsiz,mmat1,  &
       &            locmatsiz,zzero,mmat2,locmatsiz)
        
        deallocate(lbarc)

!
!       Summation over occupied bands
! 
        allocate(vect1(1:locmatsiz))
        allocate(vect2(1:locmatsiz))
        do ie1 = ibgw, nbgw
           m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
            m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
           sum=zzero
           do m = m1, m2
             do icg = 1, ncg
               ie12=(m-ibgw)*ncg+icg
               vect1(1:locmatsiz)=mmat1(1:locmatsiz,ie12)
               vect2(1:locmatsiz)=mmat2(1:locmatsiz,ie12)
               mvm=zdotc(locmatsiz,vect1,1,vect2,1)
               sum=sum+mvm
             enddo ! icg
           end do ! m
           selfxs2(ie1,ikp)=selfxs2(ie1,ikp)+sum/(m2-m1+1)
        enddo ! ie1
        deallocate(vect1,vect2)
        deallocate(mmat1,mmat2)

      endif ! wcore

!     Not needed any more ???      
      deallocate(barcs)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
      call write_cputime(fgw,tend-tstart, 'SELFXSING')
      call write_cputime(fgw,tcore-tstart,'SELFXSING: valence')
      call write_cputime(fgw,tend-tcore,  'SELFXSING: core')

      return 
      end subroutine
!EOC        
              
                 
                    
                    
                               
                        
