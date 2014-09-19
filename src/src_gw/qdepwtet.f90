!BOP
!
! !ROUTINE: qdepwtet
!
! !INTERFACE:
      subroutine qdepwtet(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the weights for q dependent BZ integrations
! using LIBBZINT
!
!
! !USES:
!
      use modmain
      use modgw
      use modmpi
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4) :: iq  !ID. number of the q-point
      
! !LOCAL VARIABLES:

      integer(4) :: ib  ! counter, run over bands
      integer(4) :: ic  ! counter, run over core states
      integer(4) :: jb  ! counter, run over bands
      integer(4) :: iom ! counter, run over frequencies
      integer(4) :: ik,ia,is,ias,sgw,ikp
      
      integer(4), allocatable :: lt(:) ! index of the q-linked tetrahedron 
     
      real(8) :: emaxb ! maximum energy of the second band
      real(8) :: edif,edsq,omsq
      
      real(8), allocatable :: bandpar(:,:)   ! Pair of bands for which the
!                                              weights are being calculated
      real(8), allocatable :: cwpar(:,:,:)   ! weight between the two bands
      real(8), allocatable :: cwparsurf(:,:,:)   ! weight between the two bands
      
      real(8) :: tstart,tend
!
!EOP
!BOC
      call cpu_time(tstart)
!
!     Initialization
!
      if (iopcore<2) unw=0.0d0
      kcw=0.0d0
      sgw=5-2*fflg
      allocate(lt(ntetnr))
      lt(1:ntetnr)=linkq(1:ntetnr,iq)
      
      if (iopcore<2) then

!---------------------------------------------------------------------
!                 core-valence     
!---------------------------------------------------------------------
      allocate(bandpar(2,nkptnr))
      allocate(cwpar(2,2,nkptnr))
      if((testid.eq.7).and.(fflg.eq.2)) then
        allocate(cwparsurf(2,2,nkptnr))
      endif
      
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
!      
!         loop over core states
!
          do ic=1,ncore(is)
!           
!           take the energies of the band
!           
            bandpar(1,1:nkptnr) = evalcr(ic,ias)
            
            do jb=numin,nstfv
!          
!             take the energies of the band
!          
              do ik=1,nkptnr
                 ikp=indkp(ik)
                 bandpar(2,ik) = evaldft(jb,ikp)
              end do
!
!             continue only if the band is at least partially unoccupied
!         
              emaxb=maxval(bandpar(2,:))
              if(emaxb.gt.efermi)then
!           
!               loop over frequencies:
!          
                do iom=1,nomeg
                  omsq=sgw*freqs(iom)*freqs(iom)
!
!                 calculate the convolution weights for the two bands
!     
                  call tetcw(nkptnr,ntetnr,2,wtetnr,bandpar,tnodesnr,lt,tvol,        &
     &                       efermi,freqs(iom),1,cwpar)
     
                  do ik=1,nkptnr
                     ikp=indkp(ik)
                     edif=evaldft(jb,ikp)-evalcr(ic,ias)
                     edsq=edif*edif
                     unw(ias,ic,jb,iom,ik)=2.0d0*cwpar(1,2,ik)*  &
     &                                     edif/(omsq-edsq)
                  enddo  
                  if((testid.eq.7).and.(fflg.eq.2))then
                    call tetcw(nkptnr,ntetnr,2,wtetnr,bandpar,tnodesnr,lt,tvol,   &
     &                       efermi,freqs(iom),4,cwpar)
                    
                    do ik=1,nkptnr
                      unwsurf(ias,ic,jb,iom,ik)=cwpar(1,2,ik)
                    enddo
                  endif
                enddo ! iom
              endif
            enddo ! jb
          enddo ! ic  
        enddo ! ia  
      enddo ! is
      deallocate(bandpar,cwpar)
      
      end if ! iopcore

!---------------------------------------------------------------------
!                 valence-valence     
!---------------------------------------------------------------------

      allocate(bandpar(nstfv,nkptnr))
      allocate(cwpar(nstfv,nstfv,nkptnr))
      if((testid.eq.7).and.(fflg.eq.2))then
        allocate(cwparsurf(nstfv,nstfv,nkptnr))
      endif
      do ik=1,nkptnr
         ikp=indkp(ik)
         bandpar(1:nstfv,ik)=evaldft(1:nstfv,ikp)
      end do

      lt(1:ntetnr)=linkq(1:ntetnr,iq)
      do iom=1,nomeg
        call tetcw(nkptnr,ntetnr,nstfv,wtetnr,bandpar,tnodesnr,lt, &
       &  tvol,efermi,freqs(iom),fflg,cwpar)
        if((testid.eq.7).and.(fflg.eq.2))then
            call tetcw(nkptnr,ntetnr,nstfv,wtetnr,bandpar,tnodesnr,lt, &
           &  tvol,efermi,freqs(iom),4,cwparsurf)
        endif
        do ik=1,nkptnr
          do ib=1,nstfv
            if(bandpar(ib,kqid(ik,iq)).gt.900.0) cwpar(1:nstfv,ib,ik)=0.0d0
            if((testid.eq.7).and.(fflg.eq.2)) then
               if(bandpar(ib,kqid(ik,iq)).gt.900.0) cwparsurf(1:nstfv,ib,ik)=0.0d0
            endif
          enddo
        enddo   
        kcw(1:nomax,numin:nstfv,iom,1:nkptnr) =           &
     &      cwpar(1:nomax,numin:nstfv,1:nkptnr)
        if((testid.eq.7).and.(fflg.eq.2)) then
          kcwsurf(1:nomax,numin:nstfv,iom,1:nkptnr) =     &
     &          cwparsurf(1:nomax,numin:nstfv,1:nkptnr)
        endif
      enddo ! iom

      deallocate(lt)
      deallocate(bandpar,cwpar)
      if((testid.eq.7).and.(fflg.eq.2)) then
        deallocate(cwparsurf)
      endif
      
      call cpu_time(tend)
      call write_cputime(fgw,tend-tstart,'QDEPWTET')

!-------------------------
!     DEBUG INFORMATION
!-------------------------
      if(debug .and. rank.eq.0)then
         if(iq==1)open(74,file='QDEPW.OUT')
         if (iopcore<2) then
         write(74,*)'------------------------------------------------------'
         write(74,*)'       convolution weights for iq =',iq
         write(74,*)'------------------------------------------------------'
         write(74,*)'------------------------------------------------------'
         write(74,*)'                   CORE '
         write(74,*)'------------------------------------------------------'
         do is=1,nspecies
           do ia=1, natoms(is)
             ias=idxas(ia,is)
             do ic=1,ncore(is)
               do ik=1,nkptnr
                 do jb=1,nstfv    
                   write(74,1)ias,ic,ik,jb,unw(ias,ic,jb,1,ik)
                 enddo      
               enddo ! jb
             enddo ! ic  
           enddo ! ia
         enddo ! is
         end if ! iopcore
         write(74,*)'------------------------------------------------------'
         write(74,*)'                   VALENCE '
         write(74,*)'------------------------------------------------------'
          do ik=1,nkptnr
           do ib=1,nstfv
             do jb=1,nstfv
               write(74,2)ik,ib,jb,kcw(ib,jb,1,ik)
             enddo
           enddo
         enddo
         if(iq==nqpt)close(74)
    1    format(' iat =',i4,' ic =',i4,' ik =',i4,' jb =',i4,' unw =',g18.10)
    2    format(' ik =',i4,' ib =',i4,' jb =',i4,' kcw =',g18.10)
      end if

      return
      end subroutine qdepwtet
!EOC          
         
