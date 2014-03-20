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

      integer(4) :: icg, ic
      integer(4) :: ie1,ie2
      integer(4) :: m,m1,m2,nmdim
      integer(4) :: ik, iq, jk
      integer(4) :: is, ia, ias

      real(8) :: tstart, tend
      real(8) :: sxs2, wkq
      
      complex(8) :: mvm     ! Sum_ij{M^i*W_ij(w)*conjg(M^j)}
      complex(8), allocatable :: minm(:,:,:) 
      complex(8), allocatable :: sxqval(:), sxqcor(:)      
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
      
      ! (non-reduced) v-diagonal basis is necessary for dealing with
      ! q->0 singularities
      call setbarcev(0.d0)
      sxs2 = -4.d0*pi*vi
      
      ! non-reduced k-grid index
      ik = idikp(ikp)
      iq = idikpq(iqp,ikp)
      jk = kqid(ik,iq)

!-------------------------------- 
!       Valence contribution
!-------------------------------- 

      allocate(minm(mbsiz,nstfv,nstfv))
      nmdim = nstfv*nstfv
      call zgemm('c','n',mbsiz,nmdim,matsiz, &
      &          zone,barcvm,matsiz,minmmat,matsiz, &
      &          zzero,minm,mbsiz)

      allocate(sxqval(ibgw:nbgw))
      sxqval(ibgw:nbgw)=zzero
     
      ! Summation over occupied bands
      do ie1 = ibgw, nbgw
         m1 = max(n12dgn(1,ie1,ikp),ibgw) ! low index
         m2 = min(n12dgn(2,ie1,ikp),nbgw) ! upper index
         sum = zzero
         do m = m1, m2
           do ie2 = 1, nomax
             mvm = zdotc(mbsiz,minm(1:mbsiz,m,ie2),1,minm(1:mbsiz,m,ie2),1)
             ! BZ integration weight
             wkq = nqptnr*wkpq(iqp,ikp)*kiw(ie2,jk)
             sum = sum+wkq*mvm
           enddo ! ie2
         end do ! m
         sxqval(ie1) = sxqval(ie1)-sum/dble(m2-m1+1)
         if (Gamma.and.(ie1.le.nomax)) then
            sxqval(ie1) = sxqval(ie1)+sxs2*singc2
         end if
      enddo ! ie1
      
      deallocate(minm)

!-------------------------------- 
!       Core electron contribution
!-------------------------------- 

      if (iopcore<=1) then
      
        allocate(minm(mbsiz,nstfv,ncg))
        nmdim = nstfv*ncg
        call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
        &          zone,barcvm,matsiz,mincmat,locmatsiz, &
        &          zzero,minm,mbsiz)

        allocate(sxqcor(ibgw:nbgw))
        sxqcor(ibgw:nbgw)=zzero

        ! Summation over occupied bands
        do ie1 = ibgw, nbgw
           sum=zzero
           m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
           m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
           do m = m1, m2
             do icg = 1, ncg
               mvm=zdotc(mbsiz,minm(1:mbsiz,m,icg),1,minm(1:mbsiz,m,icg),1)
               ! BZ integration weight
               is = corind(icg,1)
               ia = corind(icg,2)
               ias= idxas(ia,is)
               ic = corind(icg,3)
               wkq = nqptnr*wkpq(iqp,ikp)*ciw(ias,ic)
               sum = sum+wkq*mvm
             enddo ! icg
           end do ! m        
           sxqcor(ie1)=sxqcor(ie1)-sum/dble(m2-m1+1)
        enddo ! ie1
        
        deallocate(minm)
 
      endif ! core
      
!----------------------------------------
!     Sum up the contributions
!----------------------------------------

      write(96,*)'-------- CALCSELFX -------------, ikp = ', ikp, '    iqp = ', iqp
      if (iopcore.le.1) then
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
           write(96,11)ie1,sxqval(ie1),selfex(ie1,ikp)
        enddo ! ie1
        write(96,*)
        deallocate(sxqval)
      endif

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
      call write_cputime(fgw,tend-tstart, 'CALCSELFX')

   10 format(i4,2f18.10,'  | ',2f18.10,' || ',2f18.10)
   11 format(i4,2f18.10,' || ',2f18.10)
      
      return 
      end subroutine calcselfx
!EOC        
              
                 
                    
                    
                               
                        
