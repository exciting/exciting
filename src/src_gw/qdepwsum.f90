!BOP
!
! !ROUTINE: qdepw
!
! !INTERFACE:
      subroutine qdepwsum(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the weights for q dependent BZ integrations.
! \textbf{needs checking!!}
!
!
! !USES:
!
      use modmain
      use modgw
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4) :: iq  !ID. number of the q-point
      
! !LOCAL VARIABLES:

      integer(4) :: ib  ! counter, run over bands
      integer(4) :: ic  ! counter, run over core states
      integer(4) :: jb  ! counter, run over bands
      integer(4) :: iom ! counter, run over frequencies
      integer(4) :: ik,ia,is,ias,jk
      integer(4) :: ikp,jkp
      integer(4) :: fflg,sgw
      real(8) :: edif,edsq,omsq,wk
      real(8) :: tstart,tend
!
!EOP
!BOC      
      call timesec(tstart)
!
!     Initialization
!
      if (allocated(unw)) deallocate(unw) 
      allocate(unw(natmtot,ncmax,nstfv,1:freq%nomeg,nkptnr))
      unw=0.0d0

      if (allocated(kcw)) deallocate(kcw)
      allocate(kcw(nstfv,nstfv,1:freq%nomeg,nkptnr))
      kcw=0.0d0
      
      select case (freq%fconv)
        case('nofreq')
          fflg = 1
        case('refreq')
          fflg = 2
        case('imfreq')
          fflg = 3
      end select
      sgw=5-2*fflg
! 
!     loop over inequivalent atoms    
!
      wk=2.0d0/dble(kqset%nkpt)
      do ik = 1, kqset%nkpt
        jk=kqid(ik,iq)
        ikp=ik2ikp(ik)
        jkp=ik2ikp(jk)
        do ib = 1, nstfv
          if(evalsv(ib,ikp).gt.efermi)then
            if(evalsv(ib,ikp).lt.900.0)then
!------------------------------------------------------
!                   CORE-VALENCE
!------------------------------------------------------
              do is=1,nspecies
                do ia=1,natoms(is)
                  ias=idxas(ia,is)
                  do ic=1,ncore(is)
                    edif=evalsv(ib,ikp)-evalcr(ic,ias)
                    edsq=edif*edif
                    do iom=1,freq%nomeg
                      omsq=sgw*freq%freqs(iom)*freq%freqs(iom)
                      unw(ias,ic,ib,iom,ik)=cmplx(wk*edif/(omsq-edsq),0.d0)
                    enddo ! iom
                  enddo ! ic
                enddo ! ia
              enddo ! is  
            endif
          else
!------------------------------------------------------
!                   VALENCE-VALENCE
!------------------------------------------------------
            do jb=1,nstfv
              if(evalsv(jb,jkp).gt.efermi)then
                if(evalsv(jb,jkp).lt.900.0)then
                  edif=evalsv(jb,jkp)-evalsv(ib,ikp)
                  edsq=edif*edif
                  do iom=1,freq%nomeg
                    omsq=sgw*freq%freqs(iom)*freq%freqs(iom)
                    kcw(ib,jb,iom,ik)=cmplx(wk*edif/(omsq-edsq),0.d0)
                  enddo ! iom
                endif
              endif
            enddo ! jb
          endif
        enddo ! ib
      enddo ! ik  
      


!      write(74,*)'------------------------------------------------------'
!      write(74,*)'       convolution weights for iq =',iq
!      write(74,*)'------------------------------------------------------'
!      write(74,*)'------------------------------------------------------'
!      write(74,*)'                   CORE '
!      write(74,*)'------------------------------------------------------'
!      do iat=1,nat
!        do ic=1,ncore(iat)
!          do ik=1,nkpt
!            do jb=1,nb(ik)     
!              write(75,1)iat,ic,ik,jb,unw(iat,ic,jb,1,ik)
!            enddo      
!          enddo ! jb
!        enddo ! ic  
!      enddo ! iat
!      write(74,*)'------------------------------------------------------'
!      write(74,*)'                   VALENCE '
!      write(74,*)'------------------------------------------------------'
!       do ik=1,nkpt
!        do ib=1,nstfv
!          do jb=1,nstfv
!            write(75,2)ik,ib,jb,kcw(ib,jb,1,ik)
!          enddo
!        enddo
!      enddo      
!    1 format(' iat =',i4,' ic =',i4,' ik =',i4,' jb =',i4,' unw =',g18.10)
!    2 format(' ik =',i4,' ib =',i4,' jb =',i4,' kcw =',g18.10)

      call timesec(tend)
      time_bzinit = time_bzinit+tend-tstart

      end subroutine qdepwsum
!EOC          
         
