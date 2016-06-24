
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmatcorxs
! !INTERFACE:
subroutine genpmatcorxs(ik,ngp,apwalm,evecfv,evecsv,pmatc)
! !USES:
    use modinput
    use modmain
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
    integer :: i,l
    real(8) :: t1,t2
    real(8) :: kvec(3)
    real(8), dimension(nrcmtmax) :: fr,fr1,fr2,gr
    real(8), dimension(3,nrcmtmax) :: cf
!   allocatable arrays
    complex(8), allocatable :: wfmt(:,:,:)
    complex(8), allocatable :: gwfmt(:,:,:,:)
    complex(8), allocatable :: pmc(:,:,:)

!   external functions
    allocate(wfmt(lmmaxapw,nrcmtmax,nstfv))
    allocate(gwfmt(lmmaxapw,nrcmtmax,3,nstfv))
    allocate(pmc(3,ncg,nstfv))
!   set the momentum matrix elements to zero
    pmc(:,:,:)=zzero
    do i=1,3
      kvec(i)=-vkl(i,ik)
    enddo  

!   calculate momentum matrix elements in the muffin-tin
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
        do ist2=1,ncore
			l=spl(ist2,is)
            do m=1,2*spk(ist2,is)
              ist3=ist3+1
              irc=0
              do ir=1,nrmt(is),input%groundstate%lradstep
					irc=irc+1
					fr(irc)=ucore(ir,ist3)*spr(ir,is)*spr(ir,is)
              enddo ! ir
              lm1=idxlm(l,mj2ml(ist3,1))
              lm2=idxlm(l,mj2ml(ist3,2))
              do i=1,3
                do irc=1,nrcmt(is)
                  fr1(irc)=1.0d0/sqrt(2.0d0)*(preml(ist3,1)*fr(irc)*dble(gwfmt(lm1,irc,i,ist1))+ &
                  &		preml(ist3,2)*fr(irc)*dble(gwfmt(lm2,irc,i,ist1)))
                  fr2(irc)=1.0d0/sqrt(2.0d0)*(preml(ist3,1)*fr(irc)*aimag(gwfmt(lm1,irc,i,ist1))+ &
                  & 	preml(ist3,2)*fr(irc)*aimag(gwfmt(lm2,irc,i,ist1)))

!                  fr1(irc)=sqrt(3.0d0)/2.0d0*(preml(ist3,1)*fr(irc)*dble(gwfmt(lm1,irc,i,ist1))+ &
!                  &		preml(ist3,2)*fr(irc)*dble(gwfmt(lm2,irc,i,ist1)))
!                  fr2(irc)=sqrt(3.0d0)/2.0d0*(preml(ist3,1)*fr(irc)*aimag(gwfmt(lm1,irc,i,ist1))+ &
!                  & 	preml(ist3,2)*fr(irc)*aimag(gwfmt(lm2,irc,i,ist1)))

                enddo  
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
                t1=gr(nrcmt(is))
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
                t2=gr(nrcmt(is))
                pmc(i,ist3,ist1)=cmplx(t1,t2,8)
              enddo ! i
            end do ! m
          end do ! ist2
    end do ! ist1

    ! multiply by -i
    pmatc(:,:,:)=-zi*pmc(:,:,:)

    deallocate(wfmt,gwfmt,pmc)
    return
end subroutine genpmatcorxs
!EOC
