!BOP
!
! !ROUTINE: writeqp
!
! !INTERFACE:
    subroutine writeqp(signnc,znk)

! !DESCRIPTION:
! 
! This subroutine writes the qp energies to file
!
! !USES:
!
    use modmain
    use modgw
    use modmpi
      
    implicit none     

! the correlation self energy
    complex(8), intent(in) :: signnc(ibgw:nbgw,nkpt)
! delta prefactor
    real(8),    intent(in) :: znk(ibgw:nbgw,nkpt)
       
! !LOCAL VARIABLES:
      
    integer(4) :: ie   !(Counter) Runs over bands
    integer(4) :: ikp  !(Counter) Runs over k-points
      
    real(8) :: deltae,deltax
    real(8) :: ehf,eks,egw
    real(8) :: ehf0,eks0,egw0
    real(8) :: vxc,sx,sc,z

! !REVISION HISTORY:
!
! Created 16.08.05 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC

! There is a freedom in choosing the reference energy
!
    eks0 = 0.d0
    ehf0 = eks0
    egw0 = eks0
!    do ikp = 1, nkpt
!        eks0 = max(eks0,evaldft(nomax,ikp))
!        egw0 = max(egw0,eqp(nomax,ikp))
!        ehf0 = max(ehf0,evaldft(nomax,ikp)+real(selfex(nomax,ikp))-real(vxcnn(nomax,ikp)))
!    end do

!-------------------------------------------------------------------------------
    open(64,file='EVALQP.TXT',action='WRITE',form='FORMATTED')
      
      do ikp = 1, nkpt

        write(64,1) ikp, vkl(:,ikp)
        write(64,2)
        do ie = ibgw, nbgw

          eks=evaldft(ie,ikp)
          egw=eqp(ie,ikp)
          deltae=egw-eks
          
          sx=real(selfex(ie,ikp))
          sc=real(signnc(ie,ikp))
          vxc=real(vxcnn(ie,ikp))
          z=znk(ie,ikp)
          
          deltax=sx-vxc
          ehf=eks+deltax

          write(64,3) ie, eks-eks0, ehf-ehf0, egw-egw0, &
          &           sx, sc, vxc, deltax, deltae, z
        
        enddo ! ie
        write(64,*)

      enddo ! ikp
      close(64)

    1 format('k-point #',i6,':',3f12.6)
    2 format(' state    E_KS      E_HF       E_GW       Sx         Sc         Vxc         DE_HF        DE_GW       Znk')    
    3 format(i4,'  ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)

      return 
      end subroutine writeqp
!EOC      
