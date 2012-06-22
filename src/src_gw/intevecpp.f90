!BOP
!
! !ROUTINE: intevecpp
!
! !INTERFACE:
      subroutine intevecpp(ikp,jkp,ib1,ib2)

! !DESCRIPTION:
!
! This subroutine calculates the integral of the square
! modulus of the product of  two LAPW eigenvectors
!
! !USES:

      use modinput
      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ikp ! The k-point for which the (L)APW+lo
                                    ! function is ploted
      integer(4), intent(in) :: jkp ! The k-point for which the (L)APW+lo
                                    ! function is ploted
      integer(4), intent(in) :: ib1 ! The band index of the function
      integer(4), intent(in) :: ib2 ! The band index of the function
      
! !LOCAL VARIABLES:

      integer(4) :: is,ia,ias
      
      integer(4) :: ik1,ik2,ik3,ik4
      integer(4) :: iq, ir
      integer(4) :: ileb
      integer(4) :: igv(3),ig
      
      real(8) :: epii,intep

      complex(8) :: pr
      real(8) :: eprod(nrmtmax)
      real(8) :: rs(nrmtmax)
      real(8) :: gr(nrmtmax)
      real(8) :: cf(3,nrmtmax)
      real(8) :: epint(natmtot)
      
      complex(8) :: fac1,fac2,cepii
      complex(8), allocatable :: zzk(:,:), zzq(:,:)
      complex(8), allocatable :: apwalm(:,:,:,:,:)
      complex(8), allocatable :: evecmtlm1(:,:),evecmtlm2(:,:)
      complex(8), allocatable :: evecfv(:,:,:)
      complex(8), allocatable :: ylm(:)
      complex (8), Allocatable :: wfmt1(:)
      complex (8), Allocatable :: wfmt2(:)

!
! !EXTERNAL ROUTINES: 
!

! 
! !INTRINSIC ROUTINES: 
!
      intrinsic conjg
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 15th. July 2004 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC

!     Allocate the local arrays
      allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
      allocate(evecfv(nmatmax,nstfv,nspnfv))
      allocate(evecmtlm1(lmmaxapw,nrcmtmax))
      allocate(evecmtlm2(lmmaxapw,nrcmtmax))
      allocate(wfmt1(nrcmtmax))
      allocate(wfmt2(nrcmtmax))
      allocate(ylm(lmmaxapw))
      allocate(zzk(1:ngkmax,1:nstfv))
      allocate(zzq(1:ngkmax,1:nstfv))

!-----------------------------
!     MT region
!-----------------------------
      
      do is=1,nspecies

        do ia=1,natoms(is)
          ias=idxas(ia,is)

!----------------------------------------------------------------------------!
!         calculate wavefunction ( ikp, ib1 )
!----------------------------------------------------------------------------!
          call getevecfv(vklnr(:,ikp),vgklnr(:,:,:,ikp),evecfv)
          zzk(1:ngknr(1,ikp),1:nstfv)=evecfv(1:ngknr(1,ikp),1:nstfv,1)

          call match(ngknr(1,ikp),gkcnr(:,1,ikp),tpgkcnr(:,:,1,ikp), &
         &  sfacgknr(:,:,1,ikp),apwalm(:,:,:,:,1))

          call wavefmt(input%groundstate%lradstep,input%groundstate%lmaxapw, &
         &  is,ia,ngknr(1,ikp),apwalm,evecfv(:,ib1,1),lmmaxapw,evecmtlm1)

!----------------------------------------------------------------------------!
!         calculate wavefunction( jkp, ib2 )
!----------------------------------------------------------------------------!
          call getevecfv(vklnr(:,jkp),vgklnr(:,:,:,jkp),evecfv)
          zzq(1:ngknr(1,jkp),1:nstfv)=evecfv(1:ngknr(1,jkp),1:nstfv,1)
          
          call match(ngknr(1,jkp),gkcnr(:,1,jkp),tpgkcnr(:,:,1,jkp), &
         &  sfacgknr(:,:,1,jkp),apwalm(:,:,:,:,1))

          call wavefmt(input%groundstate%lradstep,input%groundstate%lmaxapw, &
         &  is,ia,ngknr(1,jkp),apwalm,evecfv(:,ib2,1),lmmaxapw,evecmtlm2)

!----------------------------------------------------------------------------!
!         calculate the product |wfmt(ikp,ib1)*conjg(wfmt(jkp,ib2))|^2
!----------------------------------------------------------------------------!
          do ir = 1, nrcmt(is)
            rs(ir)=rcmt(ir,is)*rcmt(ir,is)
          enddo
          
          eprod(:)=0.0d0
          do ileb=1,nleb
             ylm(:)=sphar(ileb,:)
          
             call zgemv('T',lmmaxapw,nrcmt(is),zone,evecmtlm1, &
            &  lmmaxapw,ylm,1,zzero,wfmt1(:),1)
             
             call zgemv('T',lmmaxapw,nrcmt(is),zone,evecmtlm2, &
            &  lmmaxapw,ylm,1,zzero,wfmt2(:),1)
          
             do ir = 1, nrcmt(is)
               pr = wfmt1(ir)*conjg(wfmt2(ir))
               eprod(ir) = eprod(ir) + &
                           4.0d0*pi*wleb(ileb)*rs(ir)*real(pr*conjg(pr))
             end do !ir

          end do !ileb
          
          call fderiv(-1,nrcmt(is),rcmt(:,is),eprod(1:nrcmt(is)),gr,cf)
          epint(ias) = gr(nrcmt(is))

          intep=intep+epint(ias)

!         output         
          write(71,11)ib1,ib2,ias,epint(ias)

        end do ! ia
      end do ! is
      deallocate(apwalm)
      deallocate(evecfv)
      deallocate(evecmtlm1)
      deallocate(evecmtlm2)
      deallocate(wfmt1)
      deallocate(wfmt2)

!-----------------------------
!     Interstitial region
!-----------------------------
      do iq=1,nqptnr
        if (kqid(ikp,iq).eq.jkp) exit
      enddo  

      cepii=zzero
      do ik1=1,ngknr(1,ikp)
        do ik2=1,ngknr(1,jkp)
          fac1=zzk(ik1,ib1)*conjg(zzq(ik2,ib2))
          do ik3=1,ngknr(1,ikp)
            fac2=fac1*conjg(zzk(ik3,ib1))
            do ik4=1,ngknr(1,jkp)
              igv(1:3)=ivg(1:3,igkig(ik1,1,ikp))- &
             &         ivg(1:3,igkig(ik2,1,jkp))- &
             &         ivg(1:3,igkig(ik3,1,ikp))+ &
             &         ivg(1:3,igkig(ik4,1,jkp))
              ig=ivgig(igv(1),igv(2),igv(3))
              cepii=cepii+fac2*zzq(ik4,ib2)*ipwint(ig)
            enddo
          enddo        
        enddo
      enddo
      epii=real(cepii)/omega

      intep=intep+epii

      write(72,10)ib1,ib2,epii
      write(73,10)ib1,ib2,intep

      deallocate(zzk)
      deallocate(zzq)

   10 format(2i4,1d15.7)
   11 format(3i4,1d15.7)
      return
    
      end subroutine intevecpp  
!EOC



