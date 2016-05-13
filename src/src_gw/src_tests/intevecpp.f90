!BOP
!
! !ROUTINE: intevecpp
!
! !INTERFACE:
      subroutine intevecpp(ik,jk,ib1,ib2)

! !DESCRIPTION:
!
! This subroutine calculates the integral of the square
! modulus of the product of  two LAPW eigenvectors
!
! !USES:

      use modinput
      use modmain
      use modgw
      use mod_ang_int

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ik  ! The k-point for which the (L)APW+lo
                                    ! function is ploted
      integer(4), intent(in) :: jk  ! The k-point for which the (L)APW+lo
                                    ! function is ploted
      integer(4), intent(in) :: ib1 ! The band index of the function
      integer(4), intent(in) :: ib2 ! The band index of the function
      
! !LOCAL VARIABLES:

      integer(4) :: is,ia,ias
      integer(4) :: ik1,ik2,ik3,ik4
      integer(4) :: iq, ir
      integer(4) :: ileb, ml, ngr
      integer(4) :: igv(3), ig
      integer(4) :: ngik, ngjk
      real(8) :: epii,intep
      complex(8) :: pr
      real(8) :: eprod(nrmtmax)
      real(8) :: rs(nrmtmax)
      real(8) :: gr(nrmtmax)
      real(8) :: cf(3,nrmtmax)
      real(8) :: epint(natmtot)
      complex(8) :: fac1,fac2,cepii
      complex(8), allocatable :: zzk(:,:), zzq(:,:)
      complex(8), allocatable :: apwalm(:,:,:,:)
      complex(8), allocatable :: evecmtlm1(:,:),evecmtlm2(:,:)
      complex(8), allocatable :: evecfv(:,:)
      complex(8), allocatable :: ylm(:)
      complex (8), Allocatable :: wfmt1(:)
      complex (8), Allocatable :: wfmt2(:)
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
      allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
      allocate(evecfv(nmatmax,nstfv))
      allocate(evecmtlm1(lmmaxapw,nrcmtmax))
      allocate(evecmtlm2(lmmaxapw,nrcmtmax))
      allocate(wfmt1(nrcmtmax))
      allocate(wfmt2(nrcmtmax))
      allocate(ylm(lmmaxapw))
      allocate(zzk(1:ngkmax,1:nstfv))
      allocate(zzq(1:ngkmax,1:nstfv))

      ngik = Gkset%ngk(1,ik)
      ngjk = Gkset%ngk(1,jk)

!---------------------------------------------
!     Initialize data for angular integration
!---------------------------------------------
      
      ml=input%groundstate%lmaxapw
      ngr=16*ml*ml/3
      call prep_ang_int(ml,ngr)
      
!-----------------------------
!     MT region
!-----------------------------
      
      do is = 1, nspecies
      do ia = 1, natoms(is)

          ias=idxas(ia,is)

!----------------------------------------------------------------------------!
!         calculate wavefunction ( ik, ib1 )
!----------------------------------------------------------------------------!
          call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik), &
          &                    nmatmax,nstsv,nspinor,evecfv)

          zzk(1:ngik,:) = evecfv(1:ngik,:)
          
          call match(ngik, &
          &          Gkset%gkc(:,1,ik), &
          &          Gkset%tpgkc(:,:,1,ik), &
          &          Gkset%sfacgk(:,:,1,ik),&
          &          apwalm)

          call wavefmt(1,input%groundstate%lmaxapw, &
          &            is,ia,ngik,apwalm,evecfv(:,ib1),lmmaxapw,evecmtlm1)

!----------------------------------------------------------------------------!
!         calculate wavefunction( jkp, ib2 )
!----------------------------------------------------------------------------!
          call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk), &
          &                    nmatmax,nstsv,nspinor,evecfv)

          zzq(1:ngjk,:) = evecfv(1:ngjk,:)
          
          call match(ngjk, &
          &          Gkset%gkc(:,1,jk), &
          &          Gkset%tpgkc(:,:,1,jk), &
          &          Gkset%sfacgk(:,:,1,jk),&
          &          apwalm)

          call wavefmt(1,input%groundstate%lmaxapw, &
          &            is,ia,ngjk,apwalm,evecfv(:,ib2),lmmaxapw,evecmtlm2)

!----------------------------------------------------------------------------!
!         calculate the product |wfmt(ikp,ib1)*conjg(wfmt(jkp,ib2))|^2
!----------------------------------------------------------------------------!
          do ir = 1, nrcmt(is)
            rs(ir) = rcmt(ir,is)*rcmt(ir,is)
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

          intep = intep+epint(ias)

!         output         
          write(71,11) ib1,ib2,ias,epint(ias)

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
      do iq = 1, kqset%nkpt
        if (kqset%kqid(ik,iq)==jk) exit
      end do

      cepii = zzero
      do ik1 = 1, ngik
        do ik2 = 1, ngjk
          fac1 = zzk(ik1,ib1)*conjg(zzq(ik2,ib2))
          do ik3 = 1, ngik
            fac2 = fac1*conjg(zzk(ik3,ib1))
            do ik4 = 1, ngjk
              igv(1:3) = ivg(1:3,Gkset%igkig(ik1,1,ik))- &
              &          ivg(1:3,Gkset%igkig(ik2,1,jk))- &
              &          ivg(1:3,Gkset%igkig(ik3,1,ik))+ &
              &          ivg(1:3,Gkset%igkig(ik4,1,jk))
              if ((igv(1).ge.intgv(1,1)).and.(igv(1).le.intgv(1,2)).and.       &
              &   (igv(2).ge.intgv(2,1)).and.(igv(2).le.intgv(2,2)).and.       &
              &   (igv(3).ge.intgv(3,1)).and.(igv(3).le.intgv(3,2)))        then
                  ig = ivgig(igv(1),igv(2),igv(3))
                  cepii = cepii+fac2*zzq(ik4,ib2)*conjg(cfunig(ig))
              end if
            enddo
          enddo        
        enddo
      enddo
      epii = real(cepii)/omega

      intep = intep+epii

      write(72,10)ib1,ib2,epii
      write(73,10)ib1,ib2,intep

      deallocate(zzk)
      deallocate(zzq)

10    format(2i4,1d15.7)
11    format(3i4,1d15.7)

      return
end subroutine intevecpp  
!EOC



