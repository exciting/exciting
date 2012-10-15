!BOP
!
! !ROUTINE: readeqp
!
! !INTERFACE:
      subroutine readeqp

! !DESCRIPTION:
! 
! This subroutine reads the qp energies from file
!
! !USES:
!
      use modmain
      use modgw
!
      implicit none  

! !LOCAL VARIABLES:

      integer(4) :: ib,nb,i,io,iu,n 
      integer(4) :: ikp          
      integer(4) :: ierr

      real(8) :: kvec(3)
      real(8) :: ehf

!      
! !REVISION HISTORY:
!
! Created 01.02.2006 by XZL
! Rewritten 2012 by DIN
!
!EOP
!BOC

      open(20,file='QPENE-eV.OUT',status='old',action='read')      
     
      read(20,*) kvec,ikp,nkp1,ib,nb
      
      if (ib>ibgw) then
        write(fgw,*)'WARNING(readepq): Mismatch ib(QPENE-eV.OUT)>ibgw.'
        write(fgw,*)'Increase ibgw from', ibgw, ' to ', ib
        ibgw=ib
      end if
      if (nb<nbgw) then
        write(fgw,*)'WARNING(readepq): Mismatch nb(QPENE-eV.OUT)<nbgw.'
        write(fgw,*)'Decrease nbgw from', nbgw, ' to ', nb
        nbgw=nb
      end if

      write(fgw,'(a,4i5)') 'readeqp: nkp1, ib, nb:', nkp1, ib, nb
      
      allocate( kvecs1(1:3,nkp1),   &
     &          eks1(ib:nb,nkp1),   &
     &          eqp1(ib:nb,nkp1),   &
     &          stat=ierr) 
      if(ierr.ne.0) then 
        write(6,*) "readeqp: fail to allocate memory"
      endif 

      rewind(20)
      
      maxoccband=0
      minunoband=nb
 
      write(fgw,'(2a5,2a16)') 'ik','ie','Enk(KS)','Enk(GW)'
      do ikp = 1, nkp1
         read(20,*) kvecs1(1:3,ikp)
         read(20,*) !skip the line
         io=0
         iu=nb
         do i = ib, nb
            read(20,*) n,eks1(i,ikp),ehf,eqp1(i,ikp)
            ! convert from eV to Ha
            eks1(i,ikp)=eks1(i,ikp)/hev
            eqp1(i,ikp)=eqp1(i,ikp)/hev
            write(fgw,'(2i5,2f16.6)') ikp,i,eks1(i,ikp),eqp1(i,ikp)
            if(eks1(i,ikp).gt.0.0d0)then
              if(i.lt.iu) iu=i
            else
              if(i.gt.io) io=i
            endif
          enddo
          if(io.gt.maxoccband) maxoccband=io
          if(iu.lt.minunoband) minunoband=iu
          write(fgw,*)
      enddo ! ikp
      close(20)

      if(maxoccband.ge.minunoband) then 
         write(fgw,'(a,3i5)') " !! maxoccband >= minunoband",maxoccband,minunoband
         write(fgw,*) "  maxoccband = minunoband - 1 is forced!"
         maxoccband = minunoband-1
      endif

      write(fgw,*)'maxoccband = ', maxoccband   
      write(fgw,*)'minunoband = ', minunoband
          
      return 
      end subroutine readeqp
!EOC      
