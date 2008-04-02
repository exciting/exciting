
! Copyright (C) 2002-2008 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpwmat
! !INTERFACE:
subroutine genpwmat(vpl,ngpmax,ngp,vgpc,gpc,igpig,ylmgp,sfacgp,vklk,ngkk, &
     igkigk,apwalmk,evecfvk,evecsvk,vklkp,ngkkp,igkigkp,apwalmkp,evecfvkp, &
     evecsvkp,pwmat)
! !USES:
  use modmain
  use modxs
!  use modxs
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pwmat  : matrix elements of plane wave (out,complex(ngq(iqcu),nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!
! !REVISION HISTORY:
!   Created November 2007 (Sagmeister)
!   Fixed bug found by Juergen Spitaler, September 2006 (JKD)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3)
  integer, intent(in) :: ngpmax
  integer, intent(in) :: ngp
  real(8), intent(in) :: vgpc(3,ngpmax)
  real(8), intent(in) :: gpc(ngpmax)
  integer, intent(in) :: igpig(ngpmax)
  complex(8), intent(in) :: ylmgp(lmmaxapw,ngpmax)
  complex(8), intent(in) :: sfacgp(ngpmax,natmtot)
  real(8), intent(in) :: vklk(3)
  integer, intent(in) :: ngkk
  integer, intent(in) :: igkigk(ngkmax)
  complex(8), intent(in) :: apwalmk(ngkmax,apwordmax,lmmaxapw,natmtot)
  complex(8), intent(in) :: evecfvk(nmatmax,nstfv)
  complex(8), intent(in) :: evecsvk(nstsv,nstsv)
  real(8), intent(in) :: vklkp(3)
  integer, intent(in) :: ngkkp
  integer, intent(in) :: igkigkp(ngkmax)
  complex(8), intent(in) :: apwalmkp(ngkmax,apwordmax,lmmaxapw,natmtot)
  complex(8), intent(in) :: evecfvkp(nmatmax,nstfv)
  complex(8), intent(in) :: evecsvkp(nstsv,nstsv)
  complex(8), intent(out) :: pwmat(ngp,nstsv,nstsv)
  ! local variables
  integer i,j,k,l,m,lm,irc,igp,ig,ir,iv(3),ivu(3),ispn,is,ia,ias,ist,jst
  integer :: igp1,igp2,ig1,ig2,iv1(3)
  real(8) :: v1(3)
  complex(8) zt1,zt2
  ! allocatable arrays
  real(8), allocatable :: jlgpr(:,:)
  complex(8), allocatable :: wfmtk(:,:,:)
  complex(8), allocatable :: wfmtkp(:,:,:)
  complex(8), allocatable :: wfmt1(:,:)
  complex(8), allocatable :: wfmt2(:,:)
  complex(8), allocatable :: zfmt(:,:)
  complex(8), allocatable :: pwfmt(:,:)
  complex(8), allocatable :: wfirk(:,:)
  complex(8), allocatable :: wfirkp(:,:)
  complex(8), allocatable :: pm(:,:,:)
  complex(8), allocatable :: cfunt(:,:), h(:,:), pmt(:,:)
  complex(8), allocatable :: evecfvt1(:,:), evecfvt2(:,:)
  ! external functions
  real(8), external :: r3taxi
  complex(8), external :: zfmtinp

  ! allocate arrays
  allocate(cfunt(ngkk,ngkkp))
  allocate(h(ngkk,nstfv))
  allocate(pmt(nstfv,nstfv))
  allocate(evecfvt1(nstfv,ngkk),evecfvt2(ngkkp,nstfv))
  ! check if q-point is commensurate with k-mesh
  if (any(abs(vpl*ngridk-nint(vpl*ngridk)).gt.epslat)) then
     write(*,*)
     write(*,'("Error(genpwmat): q-point not commensurate with k-mesh : ",&
          &3g18.10)') vpl
     write(*,*)
     stop
  end if
  ! kp-point umklapp G-vector
  v1(:)=vpl(:)+vklk(:)
  call r3frac(epslat,v1,ivu)
  ! check k+q=kp+Gw
  if (r3taxi(v1,vklkp).gt.epslat) then
     write(*,*)
     write(*,'("Error(genpwmat): q-point not commensurate with k-point and &
          &kp-point")')
     write(*,'(" k-point         : ")') vklk
     write(*,'(" q-point         : ")') vpl
     write(*,'(" kp-point        : ")') vklkp
     write(*,'(" umklapp G-vector: ")') ivu
     write(*,*)
     stop
  end if
  allocate(wfmtk(lmmaxapw,nrcmtmax,nstfv))
  allocate(wfmtkp(lmmaxapw,nrcmtmax,nstfv))
  allocate(wfmt1(lmmaxapw,nrcmtmax))
  allocate(wfmt2(lmmaxapw,nrcmtmax))
  allocate(zfmt(lmmaxapw,nrcmtmax))
  allocate(pwfmt(lmmaxapw,nrcmtmax))
  allocate(wfirk(ngrtot,nstfv))
  allocate(wfirkp(ngrtot,nstfv))
  allocate(pm(ngp,nstfv,nstfv))
  allocate(jlgpr(0:lmaxapw,nrcmtmax))
  ! zero arrays
  wfmtk(:,:,:)=zzero
  wfmtkp(:,:,:)=zzero
  wfmt1(:,:)=zzero
  wfmt2(:,:)=zzero
  zfmt(:,:)=zzero
  pwfmt(:,:)=zzero
  pm(:,:,:)=zzero
  ! loop over G+p vectors
  do igp=1,ngp
     write(*,*) 'Info(genpwmat): igp: ',igp
     ! calculate matrix elements of the plane wave in the muffin-tin
     do is=1,nspecies
        ! calculate bessel functions j_l(|G+q||r|)
        irc=0
        do ir=1,nrmt(is),lradstp
           irc=irc+1
           call sbessel(lmaxapw,gpc(igp)*spr(ir,is),jlgpr(0,irc))
           ! set up plane wave factor from Rayleigh formula
           do l=0,lmaxemat
              do m=-l,l
                 lm=idxlm(l,m)
                 pwfmt(lm,irc)=fourpi*conjg(zil(l))*jlgpr(l,irc)* &
                      conjg(ylmgp(lm,igp))
              end do
           end do
        end do
        do ia=1,natoms(is)
           ias=idxas(ia,is)
           do ist=1,nstfv
              ! calculate the wavefunction for k-point
              call wavefmt(lradstp,lmaxapw,is,ia,ngkk,apwalmk, &
                   evecfvk(1,ist),lmmaxapw,wfmtk(1,1,ist))
              ! calculate the wavefunction for kp-point
              call wavefmt(lradstp,lmaxapw,is,ia,ngkkp,apwalmkp, &
                   evecfvkp(1,ist),lmmaxapw,wfmtkp(1,1,ist))
              ! convert wavefunction and plane wave to spherical coordinates
              call zgemm('N','N',lmmaxapw,nrcmt(is),lmmaxapw,zone,zbshtapw, &
                   lmmaxapw,wfmtkp(1,1,ist),lmmaxapw,zzero,wfmt1,lmmaxapw)
              call zgemm('N','N',lmmaxapw,nrcmt(is),lmmaxapw,zone,zbshtapw, &
                   lmmaxapw,pwfmt,lmmaxapw,zzero,wfmt2,lmmaxapw)
              ! direct evaluation of exp(-i(G+q)r) on spherical grid
              ! note that the Rayleigh expansion seems to be more precise
              !do l=0,lmaxapw
              !   do m=-l,l
              !      lm=idxlm(l,m)
              !      irc=0
              !      do ir=1,nrmt(is),lradstp
              !         irc=irc+1
              !         vr(:)=spr(ir,is)*sphcov(:,lm)
              !         t1=dot_product(vgpc(:,igp),vr)
              !         wfmt2(lm,irc)=cmplx(cos(t1),-sin(t1),8)
              !      end do
              !   end do
              !end do
              ! calculate product in muffin-tin in real space
              do irc=1,nrcmt(is)
                 zfmt(:,irc)=wfmt1(:,irc)*wfmt2(:,irc)
              end do
              ! convert to spherical harmonics
              call zgemm('N','N',lmmaxapw,nrcmt(is),lmmaxapw,zone,zfshtapw, &
                   lmmaxapw,zfmt,lmmaxapw,zzero,wfmtkp(1,1,ist),lmmaxapw)
end do
           ! structure factor for G+p vector
           zt2=conjg(sfacgp(igp,ias))
           do ist=1,nstfv
              do jst=1,nstfv
                 zt1=zfmtinp(.true.,lmaxapw,nrcmt(is),rcmt(1,is),lmmaxapw, &
                      wfmtk(1,1,ist),wfmtkp(1,1,jst))
                 pm(igp,ist,jst)=pm(igp,ist,jst)+zt1*zt2
              end do
           end do
           ! end loops over atoms and species
        end do
     end do
     ! calculate matrix elements of the plane wave in the interstitial region
     forall (ist=1:nstfv)
        evecfvt1(ist,:)=conjg(evecfvk(1:ngkk,ist))
     end forall
     evecfvt2(:,:)=evecfvkp(1:ngkkp,:)
     do igp1=1,ngkk
        ig1=igkigk(igp1)
        iv1(:)=ivg(:,ig1)
        do igp2=1,ngkkp
           ig2=igkigkp(igp2)
           iv(:)=iv1(:)-ivg(:,ig2)+ivg(:,igpig(igp))+ivu(:)
           ig=ivgig(iv(1),iv(2),iv(3))
           cfunt(igp1,igp2)=cfunig(ig)
        end do
     end do
!!$h=matmul(cfunt,evecfvt2)
!!$pmt=matmul(evecfvt1,h)
     call zgemm('n','n', ngkk, nstfv, ngkkp, zone, cfunt, &
          ngkk, evecfvt2, ngkkp, zzero, h, ngkk)
     call zgemm('n','n', nstfv, nstfv, ngkk, zone, evecfvt1, &
          nstfv, h, ngkk, zzero, pmt, nstfv)
     pm(igp,:,:)=pm(igp,:,:)+pmt(:,:)
     ! compute the second-variational matrix elements of the plane wave
     if (tevecsv) then
        do i=1,nstsv
           do j=1,nstsv
              zt1=0.d0
              k=0
              do ispn=1,nspinor
                 do ist=1,nstfv
                    k=k+1
                    l=(ispn-1)*nstfv
                    do jst=1,nstfv
                       l=l+1
                       zt1=conjg(evecsvk(k,i))*evecsvkp(l,j)
                       zt1=zt1+zt1*pm(igp,ist,jst)
                    end do
                 end do
              end do
              pwmat(igp,i,j)=zt1
           end do
        end do
     else
        pwmat(:,:,:)=pm(:,:,:)
     end if
     ! end loop over G+p vectors
  end do
  deallocate(wfmtk,wfmtkp,wfmt1,wfmt2,zfmt,pwfmt,wfirk,wfirkp,pm,jlgpr)

  deallocate(cfunt,h,pmt,evecfvt1,evecfvt2)
end subroutine genpwmat
!EOC
