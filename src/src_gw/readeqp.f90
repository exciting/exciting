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
     
      read(20,*) kvec,ikp,nkp1,ibgw,nbgw
      
      allocate( kvecs1(1:3,nkp1),   &
     &          eks1(ibgw:nbgw,nkp1),   &
     &          eqp1(ibgw:nbgw,nkp1),   &
     &          stat=ierr) 
      if(ierr.ne.0) then 
        write(6,*) "readeqp: fail to allocate memory"
      endif 
 
      rewind(20)
      write(fgw,'(2a5,2a16)') 'ik','ie','Enk(KS)','Enk(GW)'
      do ikp = 1, nkp1
         read(20,*) kvecs1(1:3,ikp)
         read(20,*) !skip the line
         io=0
         iu=nb
         do i = ibgw, nbgw
            read(20,*) n,eks1(i,ikp),ehf,eqp1(i,ikp)
            ! convert from eV to Ha
            eks1(i,ikp)=eks1(i,ikp)/hev
            eqp1(i,ikp)=eqp1(i,ikp)/hev
            write(fgw,'(2i5,2f16.6)') ikp,i,eks1(i,ikp),eqp1(i,ikp)
         end do
         write(fgw,*)
      enddo ! ikp
      close(20)

      return 
      end subroutine readeqp
!EOC      
