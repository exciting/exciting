!BOP
!
! !ROUTINE: writeqp
!
! !INTERFACE:
      subroutine writeqpx

! !DESCRIPTION:
! 
! This subroutine writes the HF qp energies to file
!
! !USES:
!
      use modmain
      use modgw      

! !LOCAL VARIABLES:
      
      integer(4) :: ie   !(Counter) Runs over bands
      integer(4) :: ikp  !(Counter) Runs over k-points
      
      real(8) :: deltae,deltax
      real(8) :: ehf,eks,egw
      real(8) :: vxc,sx,sc,z
      real(8), dimension(3) :: kvec
      real(8) :: evmks,evmhf,evmgw

! !REVISION HISTORY:
!
! Created 16.08.05 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC
      evmks=maxval(evaldft(maxoccband,:))
      evmhf=maxval(evaldft(maxoccband,:)+real(selfex(maxoccband,:))-real(vxcnn(maxoccband,:)))
      evmgw=maxval(eqp(maxoccband,:))

      open(64,file='QPENE-eV.OUT',action='WRITE',form='FORMATTED')
      do ikp = 1, nkpt
        kvec(1:3)=vkl(1:3,ikp)
        write(64,1)kvec,ikp,nkpt,ibgw,nbgw
        write(64,2)
        do ie = ibgw, nbgw
          
          eks=evaldft(ie,ikp)-evmks
          
          deltax=real(selfex(ie,ikp))-real(vxcnn(ie,ikp))
          ehf=evaldft(ie,ikp)+deltax-evmhf
          
          egw=eqp(ie,ikp)-evmgw         

          deltae=eqp(ie,ikp)-evaldft(ie,ikp)
          
          sx=real(selfex(ie,ikp))
          sc=0.d0
          vxc=real(vxcnn(ie,ikp))
          
          z=0.d0

          write(64,3) ie,eks*hev,ehf*hev,egw*hev, &
         &            sx*hev,sc*hev,vxc*hev,      &
         &            deltax*hev,deltae*hev,z
        enddo
        write(64,*)
      enddo    
      close(64)    
    1 format(3e19.12,4i6)
    2 format(' state    E_KS      E_HF       E_GW       Sx         Sc         Vxc         DE_HF        DE_GW       Znk')    
    3 format(i4,'  ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)     

      return 
      end subroutine writeqpx
!EOC      
