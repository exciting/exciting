
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmat
! !INTERFACE:
subroutine genpmat(ngp,igpig,vgpc,apwalm,evecfv,evecsv,pmat)
! !USES:
use modmain
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
!   Created November 2003 (Sharma)
!   Fixed bug found by Juergen Spitaler, September 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: pmat(3,nstsv,nstsv)
! local variables
integer ispn,is,ia,ist,jst
integer i,j,k,l,igp,ifg,ir
complex(8) zsum,zt1,zv(3)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: gwfmt(:,:,:,:)
complex(8), allocatable :: wfir(:,:)
complex(8), allocatable :: gwfir(:,:,:)
complex(8), allocatable :: pm(:,:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(wfmt(lmmaxapw,nrcmtmax,nstfv))
allocate(gwfmt(lmmaxapw,nrcmtmax,3,nstfv))
allocate(wfir(ngrtot,nstfv))
allocate(gwfir(ngrtot,3,nstfv))
allocate(pm(3,nstfv,nstfv))
! set the momentum matrix elements to zero
pm(:,:,:)=0.d0
! calculate momentum matrix elements in the muffin-tin
do is=1,nspecies
  do ia=1,natoms(is)
    do ist=1,nstfv
! calculate the wavefunction
      call wavefmt(lradstp,lmaxapw,is,ia,ngp,apwalm,evecfv(:,ist),lmmaxapw, &
       wfmt(:,:,ist))
! calculate the gradient
      call gradzfmt(lmaxapw,nrcmt(is),rcmt(:,is),lmmaxapw,nrcmtmax, &
       wfmt(:,:,ist),gwfmt(:,:,:,ist))
    end do
    do ist=1,nstfv
      do jst=ist,nstfv
        do i=1,3
          zt1=zfmtinp(.true.,lmaxapw,nrcmt(is),rcmt(:,is),lmmaxapw, &
           wfmt(:,:,ist),gwfmt(:,:,i,jst))
          pm(i,ist,jst)=pm(i,ist,jst)+zt1
        end do
      end do
    end do
  end do
end do
! calculate momemntum matrix elements in the interstitial region
wfir(:,:)=0.d0
gwfir(:,:,:)=0.d0
do ist=1,nstfv
  do igp=1,ngp
    ifg=igfft(igpig(igp))
    zt1=evecfv(igp,ist)
    wfir(ifg,ist)=zt1
! calculate the gradient
    do i=1,3
      gwfir(ifg,i,ist)=zi*vgpc(i,igp)*zt1
    end do
  end do
! Fourier transform the wavefunction to real-space
  call zfftifc(3,ngrid,1,wfir(:,ist))
  do i=1,3
    call zfftifc(3,ngrid,1,gwfir(:,i,ist))
  end do
end do
! find the overlaps
do ist=1,nstfv
  do jst=ist,nstfv
    do i=1,3
      zsum=0.d0
      do ir=1,ngrtot
        zsum=zsum+cfunir(ir)*conjg(wfir(ir,ist))*gwfir(ir,i,jst)
      end do
      zt1=zsum/dble(ngrtot)
      pm(i,ist,jst)=pm(i,ist,jst)+zt1
    end do
  end do
end do
! multiply by -i and set lower triangular part
do ist=1,nstfv
  do jst=ist,nstfv
    pm(:,ist,jst)=-zi*pm(:,ist,jst)
    pm(:,jst,ist)=conjg(pm(:,ist,jst))
  end do
end do
! compute the second-variational momentum matrix elements
if (tevecsv) then
  do i=1,nstsv
    do j=1,nstsv
      zv(:)=0.d0
      k=0
      do ispn=1,nspinor
        do ist=1,nstfv
          k=k+1
          l=(ispn-1)*nstfv
          do jst=1,nstfv
            l=l+1
            zt1=conjg(evecsv(k,i))*evecsv(l,j)
            zv(:)=zv(:)+zt1*pm(:,ist,jst)
          end do
        end do
      end do
      pmat(:,i,j)=zv(:)
    end do
  end do
else
  pmat(:,:,:)=pm(:,:,:)
end if
deallocate(wfmt,gwfmt,wfir,gwfir,pm)
return
end subroutine
!EOC

