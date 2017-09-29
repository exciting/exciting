
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmatcorxs
! !INTERFACE:
subroutine genpmatcorxs(ik,ngp,apwalm,evecfv,evecsv,pmatc)
! !USES:
    use modinput
    use mod_constants, only: zzero, zi
    use mod_kpoint, only: vkl
    use mod_atoms, only: idxas, natmtot
    use mod_muffin_tin, only: nrcmtmax, nrcmt, rcmt, nrmtmax, nrmt, &
      & idxlm, lmmaxapw
    use mod_atoms, only: spl, spk, spr
    use mod_spin, only: nspinor
    use mod_Gkvector, only: ngkmax
    use mod_APW_LO, only: apwordmax
    use mod_eigensystem, only: nmatmax 
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use m_b_getgrst, only: b_wavefmtsv
    !use modmain
    use modxas
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point position (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmatc   : momentum matrix elements (out,complex(3,ncg,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements between a core state and a conduction
!   state.
!   
!
! !REVISION HISTORY:
!   Created August 2006 (RGA)
!   Revisited June 2011 (DIN)
!	Adjusted to core states November 2015 (Christian Vorwerk)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in)     :: ik
    integer, intent(in)     :: ngp
    complex(8), intent(in)  :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: pmatc(3,ncg,nstsv)
!   local variables
    integer :: is,ia,ist1,ist2,ist3,ias,lm1,lm2,m,ir,irc
    integer :: i,l, ispn
    real(8) :: t1,t2
    real(8) :: kvec(3)
    real(8), dimension(nrcmtmax) :: fr,fr1,fr2,gr
    real(8), dimension(3,nrcmtmax) :: cf
!   allocatable arrays
    complex(8), allocatable :: wfmt(:,:,:)
    complex(8), allocatable :: gwfmt(:,:,:,:)
    complex(8), allocatable :: pmc(:,:,:)
    complex(8), allocatable :: wfmtsv(:,:,:), gwfmtsv(:,:,:,:,:)

    if (.not. input%groundstate%tevecsv) then
!     external functions
      allocate(wfmt(lmmaxapw,nrcmtmax,nstfv))
      allocate(gwfmt(lmmaxapw,nrcmtmax,3,nstfv))
      allocate(pmc(3,ncg,nstfv))
!     set the momentum matrix elements to zero
      pmc(:,:,:)=zzero
      do i=1,3
        kvec(i)=-vkl(i,ik)
      enddo  
!     calculate momentum matrix elements in the muffin-tin
      do ist1=1,nstfv 
        ist3=0
        ia=input%xs%bse%xasatom
        is=input%xs%bse%xasspecies
        ias=idxas(ia,is)
        ! calculate the wavefunction
        call wavefmt(input%groundstate%lradstep, &
          &  input%groundstate%lmaxapw,is,ia,ngp,apwalm, &
          &  evecfv(:,ist1),lmmaxapw,wfmt(:,:,ist1))
        ! calculate the gradient
        call gradzfmt(input%groundstate%lmaxapw,nrcmt(is), &
          &  rcmt(:,is),lmmaxapw,nrcmtmax,wfmt(:,:,ist1),    &
          &  gwfmt(:,:,:,ist1))
        do ist2=xasstart,xasstop
          irc=0
          do ir=1,nrmt(is),input%groundstate%lradstep
            irc=irc+1
            fr(irc)=ucore(ir,ist2)*spr(ir,is)*spr(ir,is)
          enddo ! ir
          lm1=idxlm(lxas,mj2ml(lxas,mj(ist2),1))
          lm2=idxlm(lxas,mj2ml(lxas,mj(ist2),2))
          do i=1,3
            do irc=1,nrcmt(is)
              fr1(irc)=1.0d0/sqrt(2.0d0)*(preml(lxas,spj(ist2),mj(ist2),1)*fr(irc)& 
                &   *dble(gwfmt(lm1,irc,i,ist1))+ preml(lxas,spj(ist2),mj(ist2),2)&
                &   *fr(irc)*dble(gwfmt(lm2,irc,i,ist1)))
               fr2(irc)=1.0d0/sqrt(2.0d0)*(preml(lxas,spj(ist2),mj(ist2),1)*fr(irc)&
                 &   *aimag(gwfmt(lm1,irc,i,ist1))+ preml(lxas,spj(ist2),mj(ist2),2)&
                 &   *fr(irc)*aimag(gwfmt(lm2,irc,i,ist1)))
            enddo  
            call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
            t1=gr(nrcmt(is))
            call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
            t2=gr(nrcmt(is))
            pmc(i,ist2,ist1)=cmplx(t1,t2,8)
          enddo ! i
        end do ! ist2
      end do ! ist1

      ! multiply by -i
      pmatc(:,:,:)=-zi*pmc(:,:,:)
      deallocate(wfmt,gwfmt,pmc)
    else
!     external functions
      allocate(gwfmtsv(lmmaxapw,nrmtmax,3,nspinor,nstsv))
      allocate(wfmtsv(lmmaxapw,nrmtmax,nspinor))
      allocate(pmc(3,ncg,nstsv))
      do ist1=1,nstsv 
        ist3=0
        ia=input%xs%bse%xasatom
        is=input%xs%bse%xasspecies
        ias=idxas(ia,is)
        ! calculate the gradient
        call b_wavefmtsv(input%groundstate%lradstep,input%groundstate%lmaxapw,is ,ia , ngp,&
         &   ist1, apwalm, evecfv, evecsv, wfmtsv)
        do ispn=1,nspinor
          call gradzfmt(input%groundstate%lmaxapw,nrcmt(is), &
          &  rcmt(:,is),lmmaxapw,nrmtmax,wfmtsv(:,:,ispn),    &
          &  gwfmtsv(:,:,:,ispn,ist1))
        end do
        do ist2=xasstart, xasstop
          irc=0
          do ir=1,nrmt(is),input%groundstate%lradstep
            irc=irc+1
            fr(irc)=ucore(ir,ist2)*spr(ir,is)*spr(ir,is)
          enddo ! ir
          lm1=idxlm(lxas,mj2ml(lxas,mj(ist2),1))
          lm2=idxlm(lxas,mj2ml(lxas,mj(ist2),2))
          do i=1,3
            do irc=1,nrmt(is)
              fr1(irc)=preml(lxas,spj(ist2),mj(ist2),1)*fr(irc)*dble(gwfmtsv(lm1,irc,i,1,ist1))+ &
                &   preml(lxas,spj(ist2),mj(ist2),2)*fr(irc)*dble(gwfmtsv(lm2,irc,i,2,ist1))
              fr2(irc)=preml(lxas,spj(ist2),mj(ist2),1)*fr(irc)*aimag(gwfmtsv(lm1,irc,i,1,ist1))+ &
                &   preml(lxas,spj(ist2),mj(ist2),2)*fr(irc)*aimag(gwfmtsv(lm2,irc,i,2,ist1))
            enddo  
            call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
            t1=gr(nrcmt(is))
            call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
            t2=gr(nrcmt(is))
            pmc(i,ist2,ist1)=cmplx(t1,t2,8)
          enddo ! i
        end do ! ist2
      end do ! ist1

      ! multiply by -i
      pmatc(:,:,:)=-zi*pmc(:,:,:)
      deallocate(gwfmtsv,pmc, wfmtsv)
    end if
    return
end subroutine genpmatcorxs
!EOC
