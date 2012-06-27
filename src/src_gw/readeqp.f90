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

      integer(4) :: ib,i,io,iu   !(Counter) Runs over bands
      integer(4) :: ikp          !(Counter) Runs over k-points
      integer(4) :: ierr

      real(8) :: kvec(3)
      real(8) :: ehf

!      
! !REVISION HISTORY:
!
! Created 01.02.2006 by XZL
!
!EOP
!BOC

      open(20,file='QPENE-eV.OUT',status='old',action='read')      
     
      read(20,*) kvec,ikp,nkp1,ibgw,nbgw
      nbandsgw=nbgw-ibgw+1

      write(fgw,'(a,4i5)') 'readeqp: nkp1, ibgw, nbgw, nbandsgw :', nkp1, ibgw, nbgw, nbandsgw
      
      allocate( kvecs1(1:3,nkp1),        &
     &          eks1(ibgw:nbgw,nkp1),   &
     &          eqp1(ibgw:nbgw,nkp1),   &
     &          stat=ierr) 
      if(ierr.ne.0) then 
        write(6,*) "readeqp: fail to allocate memory"
      endif 

      rewind(20)
      
      maxoccband=0
      minunoband=nbgw
 
      write(fgw,'(2a5,2a16)') 'ik','ie','Enk(KS)','Enk(GW)'
      do ikp = 1, nkp1
         read(20,*) kvecs1(1:3,ikp)
         read(20,*) !skip the line
         io=0
         iu=nbgw
         do ib = ibgw, nbgw
            read(20,*) i,eks1(ib,ikp),ehf,eqp1(ib,ikp)
            ! convert from eV to Ha
            eks1(ib,ikp)=eks1(ib,ikp)/hev
            eqp1(ib,ikp)=eqp1(ib,ikp)/hev
            write(fgw,'(2i5,2f16.6)') ikp,ib,eks1(ib,ikp),eqp1(ib,ikp)
            if(eks1(ib,ikp).gt.0.0d0)then
              if(ib.lt.iu) iu=ib
            else
              if(ib.gt.io) io=ib
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
