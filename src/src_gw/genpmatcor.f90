
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmat
! !INTERFACE:
subroutine genpmatcor(ik,ngp,apwalm,evecfv,evecsv,pmatc)
! !USES:
    use modinput
    use modmain
    use modgw
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmat   : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!
! !REVISION HISTORY:
!   Created August 2006 (RGA)
!   Revisited June 2011 (DIN)
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
    integer :: is,ia,ist1,ist2,ist3,ias,lm,m,ir,irc
    integer :: i,l
    real(8) :: arg,t1,t2
    real(8) :: kvec(3)
    real(8), dimension(nrcmtmax) :: fr,fr1,fr2,gr
    real(8), dimension(3,nrcmtmax) :: cf
    complex(8) :: phs
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
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          arg=atposl(1,ia,is)*kvec(1)+ &
         &    atposl(2,ia,is)*kvec(2)+ &
         &    atposl(3,ia,is)*kvec(3)
          phs=cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)
          ! calculate the wavefunction
          call wavefmt(input%groundstate%lradstep, &
         &  input%groundstate%lmaxapw,is,ia,ngp,apwalm, &
         &  evecfv(:,ist1),lmmaxapw,wfmt(:,:,ist1))
          ! calculate the gradient
          call gradzfmt(input%groundstate%lmaxapw,nrcmt(is), &
         &  rcmt(:,is),lmmaxapw,nrcmtmax,wfmt(:,:,ist1),    &
         &  gwfmt(:,:,:,ist1))
          do ist2=1,ncore(is)
            l=spl(ist2,is)
            irc=0
            do ir=1,nrmt(is),input%groundstate%lradstep
                irc=irc+1
                fr(irc)=ucore(ir,1,ist2,ias)*spr(ir,is)*spr(ir,is)
            enddo ! ir
            do m=-l,l
              lm=idxlm(l,m)
              ist3=ist3+1
              do i=1,3
                do irc=1,nrcmt(is)
                  fr1(irc)=fr(irc)*dble(gwfmt(lm,irc,i,ist1))
                  fr2(irc)=fr(irc)*aimag(gwfmt(lm,irc,i,ist1))
                enddo  
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
                t1=gr(nrcmt(is))
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
                t2=gr(nrcmt(is))
                pmc(i,ist3,ist1)=cmplx(t1,t2,8)*phs
              enddo ! i
            end do ! m
          end do ! ist2
        end do ! ia
      end do ! is
    end do ! ist1

    ! multiply by -i
    pmatc(:,:,:)=-zi*pmc(:,:,:)

    deallocate(wfmt,gwfmt,pmc)
    return
end subroutine
!EOC
