!BOP
!
! !ROUTINE: selfesing
!
! !INTERFACE: 
      subroutine selfcsing(ikp)
      
! !DESCRIPTION:
!
! This subroutine calculates the singular terms of the selfenergy.
!
! !USES:

      use modmain
      use modgw

! !LOCAL VARIABLES:            

      implicit none
      integer(4), intent(in) :: ikp

      integer(4) :: ia 
      integer(4) :: is
      integer(4) :: ic      ! (Counter) Runs over core states.
      integer(4) :: ias     ! (Counter) Runs over all atoms
      integer(4) :: ie1     ! (Counter) Runs over bands.
      integer(4) :: ie12
      integer(4) :: ie2     ! (Counter) Runs over bands.
      integer(4) :: iom     ! (Counter) Runs over frequencies.
      integer(4) :: m, m1, m2
      integer(4) :: dimtk
      
      real(8)    :: tstart, tend, tcore
      real(8)    :: enk
      
      complex(8) :: sum1, sum2
      complex(8), allocatable :: mw1m(:,:,:)     ! Sum_ij{M^i*Ws1_ij*conjg(M^j)}
      complex(8), allocatable :: mw2m(:,:,:)     ! Sum_ij{M^i*Ws2_ij*conjg(M^j)}
      complex(8), allocatable :: mmat(:,:)
      complex(8), allocatable :: wmat(:,:), tmat1(:,:), tmat2(:,:)
      complex(8), allocatable :: sc1(:), sc2(:)

! !INTRINSIC ROUTINES: 

      
      intrinsic abs
      intrinsic sign
      intrinsic cmplx
      intrinsic conjg
      intrinsic cpu_time
      

! !EXTERNAL ROUTINES: 
      
      complex(8), external :: zdotc
      external zgemv
      complex(8), external :: freqconvl
 
! !REVISION HISTORY:
!
! Created 23.06.05 by RGA.
! Revisted July 2011 by DIN
!
!EOP
!BOC
!
      call cpu_time(tstart)
      
!-------------------------------- 
!       Valence contribution
!-------------------------------- 

      dimtk=nbandsgw*nstfv
      allocate(mmat(1:matsiz,1:dimtk))
      ie12=0
      do ie1 = ibgw, nbgw
        do ie2 = 1, nstfv
           ie12=ie12+1
           mmat(1:matsiz,ie12)=minmmat(1:matsiz,ie1,ie2)
        end do
      end do  

      allocate(mw1m(nomeg,ibgw:nbgw,nstfv))
      allocate(mw2m(nomeg,ibgw:nbgw,nstfv))
      allocate(wmat(matsiz,matsiz))
      allocate(tmat1(matsiz,dimtk))
      allocate(tmat2(matsiz,dimtk))
      do iom = 1, nomeg
        
        wmat=ws1(1:matsiz,1:matsiz,iom)
        call zhemm('l','u',matsiz,dimtk, &
       &  zone,wmat,matsiz,mmat,matsiz,zzero,tmat1,matsiz)
        
        wmat=ws2(1:matsiz,1:matsiz,iom)
        call zhemm('l','u',matsiz,dimtk, &
       &  zone,wmat,matsiz,mmat,matsiz,zzero,tmat2,matsiz)
       
        ie12=0
        do ie1 = ibgw, nbgw
          do ie2 = 1, nstfv
            ie12=ie12+1
            mw1m(iom,ie1,ie2)=zdotc(matsiz,mmat(:,ie12),1,tmat1(:,ie12),1)
            mw2m(iom,ie1,ie2)=zdotc(matsiz,mmat(:,ie12),1,tmat2(:,ie12),1)
          end do
        end do
      end do !iom
      deallocate(mmat,wmat,tmat1,tmat2)
      
      allocate(sc1(ibgw:nbgw))
      allocate(sc2(ibgw:nbgw))
      do iom = 1, nomeg
        ! Sum over ie2
        sc1(:)=zzero
        sc2(:)=zzero
        do ie1 = ibgw, nbgw
          do ie2 = 1, nstfv
            enk=evaldft(ie2,ikp)
            sc1(ie1)=sc1(ie1)+freqconvl(iom,nomeg,freqs(iom),enk,mw1m(1:nomeg,ie1,ie2),freqs,womeg)
            sc2(ie1)=sc2(ie1)+freqconvl(iom,nomeg,freqs(iom),enk,mw2m(1:nomeg,ie1,ie2),freqs,womeg)
          end do ! ie2
        enddo ! ie1
        ! Average over the degenerate states
        do ie1 = ibgw, nbgw
          m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
          m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
          sum1=zzero
          sum2=zzero
          do m = m1, m2
            sum1=sum1+sc1(m)
            sum2=sum2+sc2(m)
          end do ! m
          selfcs1(ie1,ikp,iom)=selfcs1(ie1,ikp,iom)+sum1/(m2-m1+1)
          selfcs2(ie1,ikp,iom)=selfcs2(ie1,ikp,iom)+sum2/(m2-m1+1)
        end do ! ie1
      end do ! iom
      deallocate(sc1,sc2)
      deallocate(mw1m,mw2m)

!------------------------------------
!     Core-valence contribution
!------------------------------------

      call cpu_time(tcore)

      if(wcore)then      

        dimtk=nbandsgw*ncg
        allocate(mmat(locmatsiz,dimtk))
        ie12=0
        do ie1 = ibgw, nbgw
          do ie2 = 1, ncg
             ie12=ie12+1
             mmat(1:locmatsiz,ie12)=mincmat(1:locmatsiz,ie1,ie2)
          end do
        end do  
!
!       Sum_ij{M^i*W^c_{ij}*conjg(M^j)}
!
        allocate(mw1m(nomeg,ibgw:nbgw,ncg))
        allocate(mw2m(nomeg,ibgw:nbgw,ncg))
        allocate(wmat(locmatsiz,locmatsiz))
        allocate(tmat1(locmatsiz,dimtk))
        allocate(tmat2(locmatsiz,dimtk))
        do iom = 1, nomeg

          wmat=ws1(1:locmatsiz,1:locmatsiz,iom)
          call zhemm('l','u',locmatsiz,dimtk, &
         &  zone,wmat,locmatsiz,mmat,locmatsiz,zzero,tmat1,locmatsiz)
        
          wmat=ws2(1:locmatsiz,1:locmatsiz,iom)
          call zhemm('l','u',locmatsiz,dimtk, &
         &  zone,wmat,locmatsiz,mmat,locmatsiz,zzero,tmat2,locmatsiz)
         
          ie12=0
          do ie1 = ibgw, nbgw
            do ie2 = 1, ncg
              ie12=ie12+1
              mw1m(iom,ie1,ie2)=zdotc(locmatsiz,mmat(:,ie12),1,tmat1(:,ie12),1)
              mw2m(iom,ie1,ie2)=zdotc(locmatsiz,mmat(:,ie12),1,tmat2(:,ie12),1)
            end do
          end do
        end do !iom
        deallocate(mmat,wmat,tmat1,tmat2)

        allocate(sc1(ibgw:nbgw))
        allocate(sc2(ibgw:nbgw))
        do iom = 1, nomeg
          ! Sum over ie2
          sc1(:)=zzero
          sc2(:)=zzero
          do ie1 = ibgw, nbgw
            do ie2 = 1, ncg
              is=corind(ie2,1)
              ia=corind(ie2,2)
              ias=idxas(ia,is)
              ic=corind(ie2,3)
              enk=evalcr(ic,ias)
              sc1(ie1)=sc1(ie1)+freqconvl(iom,nomeg,freqs(iom),enk,mw1m(1:nomeg,ie1,ie2),freqs,womeg)
              sc2(ie1)=sc2(ie1)+freqconvl(iom,nomeg,freqs(iom),enk,mw2m(1:nomeg,ie1,ie2),freqs,womeg)
            end do ! ie2
          enddo ! ie1
          ! Average over the degenerate states
          do ie1 = ibgw, nbgw
            m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
            m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
            sum1=zzero
            sum2=zzero
            do m = m1, m2
              sum1=sum1+sc1(m)
              sum2=sum2+sc2(m)
            end do ! m
            selfcs1(ie1,ikp,iom)=selfcs1(ie1,ikp,iom)+sum1/(m2-m1+1)
            selfcs2(ie1,ikp,iom)=selfcs2(ie1,ikp,iom)+sum2/(m2-m1+1)
          end do ! ie1
        enddo ! iom
        deallocate(sc1,sc2)
        deallocate(mw1m,mw2m)

      endif ! wcore

!----------------------------------------
!     Deallocate the global array
!----------------------------------------
      deallocate(ws1,ws2)   

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

      call write_cputime(fgw,tend-tstart, 'SELFCSING')
      call write_cputime(fgw,tcore-tstart,'SELFCSING: valence')
      call write_cputime(fgw,tend-tcore,  'SELFCSING: core')

      return 
      end subroutine selfcsing
!EOC        
              
                 
                    
                    
                               
                        
