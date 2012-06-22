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
      integer(4) :: ik,ia,is,ias,sgw,jk
      integer(4) :: ikp,jkp
      
      real(8) :: edif,edsq,omsq,wk
      real(8) :: tstart,tend
!
!EOP
!BOC      
      call cpu_time(tstart)
! 
!     loop over inequivalent atoms    
!
      sgw=5-2*fflg
      wk=2.0d0/dble(nqptnr)
      do ik = 1, nqptnr
        jk=kqid(ik,iq)
        ikp=indkp(ik)
        jkp=indkp(jk)
        do ib = 1, nstfv
          if(evaldft(ib,ikp).gt.efermi)then
            if(evaldft(ib,ikp).lt.900.0)then
!------------------------------------------------------
!                   CORE-VALENCE
!------------------------------------------------------
              do is=1,nspecies
                do ia=1,natoms(is)
                  ias=idxas(ia,is)
                  do ic=1,ncore(is)
                    edif=evaldft(ib,ikp)-evalcr(ic,ias)
                    edsq=edif*edif
                    do iom=1,nomeg
                      omsq=sgw*freqs(iom)*freqs(iom)
                      unw(ias,ic,ib,iom,ik)=wk*edif/(omsq-edsq)
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
              if(evaldft(jb,jkp).gt.efermi)then
                if(evaldft(jb,jkp).lt.900.0)then
                  edif=evaldft(jb,jkp)-evaldft(ib,ikp)
                  edsq=edif*edif
                  do iom=1,nomeg
                    omsq=sgw*freqs(iom)*freqs(iom)
                    kcw(ib,jb,iom,ik)=wk*edif/(omsq-edsq)
                  enddo ! iom
                endif
              endif
            enddo ! jb
          endif
        enddo ! ib
      enddo ! ik  
      
      call cpu_time(tend)
      call write_cputime(fgw,tend-tstart,'QDEPWSUM')


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

      end subroutine qdepwsum
!EOC          
         
