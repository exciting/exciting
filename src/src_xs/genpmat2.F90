

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmat2
! !INTERFACE:


subroutine genpmat2(ngp, igpig, vgpc, evecfv, evecsv, pmat)
! !USES:
use modinput
  use modmain
  use modxs, only: apwcmt, locmt, ripaa, ripalo, riploa, riplolo
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmat   : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!   The gradient is applied explicitly only to the radial functions and 
!   corresponding spherical harmonics for the muffin-tin part. In the 
!   interstitial region the gradient is evaluated analytically.
!   Parts taken from the routine {\tt genpmat}.
!
! !REVISION HISTORY:
!   Created April 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: ngp
  integer, intent(in) :: igpig(ngkmax)
  real(8), intent(in) :: vgpc(3, ngkmax)
  complex(8), intent(in) :: evecfv(nmatmax, nstfv)
  complex(8), intent(in) :: evecsv(nstsv, nstsv)
  complex(8), intent(out) :: pmat(3, nstsv, nstsv)
  ! local variables
  integer :: ispn, is, ia, ias, ist, jst
  integer :: ist1, l1, m1, lm1, l3, m3, lm3, io, io1, io2, ilo, ilo1, ilo2
  integer :: i, j, k, l
  integer :: igp1, igp2, ig1, ig2, ig, iv1(3), iv(3)
  complex(8) :: zt1, zv(3)
  ! allocatable arrays
  complex(8), allocatable :: wfmt(:, :, :)
  complex(8), allocatable :: gwfmt(:, :, :, :)
  complex(8), allocatable :: pm(:, :, :)
  complex(8), allocatable :: cfunt(:, :), h(:, :), pmt(:, :)
  complex(8), allocatable :: evecfv1(:, :), evecfv2(:, :)
  complex(8), allocatable :: zv2(:)
  ! external functions
  complex(8) zfmtinp
  external zfmtinp
  allocate(zv2(nstfv))
  allocate(wfmt(lmmaxapw, nrcmtmax, nstfv))
  allocate(gwfmt(lmmaxapw, nrcmtmax, 3, nstfv))
  allocate(cfunt(ngp, ngp))
  allocate(h(ngp, nstfv))
  allocate(pmt(nstfv, nstfv))
  allocate(evecfv1(nstfv, ngp), evecfv2(ngp, nstfv))
  allocate(pm(nstfv, nstfv, 3))
  ! set the momentum matrix elements to zero
  pm(:, :, :)=0.d0
  ! loop over species and atoms
  do is=1, nspecies
     do ia=1, natoms(is)
	ias=idxas(ia, is)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
	do j=1, 3
	   do l1=0, input%groundstate%lmaxapw
	      do m1=-l1, l1
		 lm1=idxlm(l1, m1)
		 do io1=1, apword(l1, is)
		    zv2(:)=zzero
		    do l3=0, input%groundstate%lmaxapw
		       do m3=-l3, l3
			  lm3=idxlm(l3, m3)
			  do io2=1, apword(l3, is)
			     call zaxpy(nstfv, &
				  zone * ripaa(io1, lm1, io2, lm3, ias, j), &
				  apwcmt(1, io2, lm3, ias), 1, zv2, 1)
			  end do
		       end do
		    end do
		    call zoutpr(nstfv, nstfv, &
			 zone, apwcmt(1, io1, lm1, ias), zv2, pm(1, 1, j))
		 end do
	      end do
	   end do
	end do
	if (nlotot.gt.0) then
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
	   do j=1, 3
	      do l3=0, input%groundstate%lmaxapw
		 do m3=-l3, l3
		    lm3=idxlm(l3, m3)
		    do io=1, apword(l3, is)
		       zv2(:)=zzero
		       do ilo=1, nlorb(is)
			  l1=lorbl(ilo, is)
			  do m1=-l1, l1
			     lm1=idxlm(l1, m1)
			     call zaxpy(nstfv, &
				  zone * ripalo(io, lm3, ilo, m1, ias, j), &
				  locmt(1, ilo, m1, ias), 1, zv2, 1)
			  end do
		       end do
		       call zoutpr(nstfv, nstfv, &
			    zone, apwcmt(1, io, lm3, ias), zv2, pm(1, 1, j))
		    end do
		 end do
	      end do
	   end do
           !--------------------------------------!
           !     local-orbital-APW contribution   !
           !--------------------------------------!
	   do j=1, 3
	      do ilo=1, nlorb(is)
		 l1=lorbl(ilo, is)
		 do m1=-l1, l1
		    lm1=idxlm(l1, m1)
		    zv2(:)=zzero
		    do l3=0, input%groundstate%lmaxapw
		       do m3=-l3, l3
			  lm3=idxlm(l3, m3)
			  do io=1, apword(l3, is)
			     call zaxpy(nstfv, &
				  zone * riploa(ilo, m1, io, lm3, ias, j), &
				  apwcmt(1, io, lm3, ias), 1, zv2, 1)
			  end do
		       end do
		    end do
		    call zoutpr(nstfv, nstfv, &
			 zone, locmt(1, ilo, m1, ias), zv2, pm(1, 1, j))
		 end do
	      end do
	   end do
           !------------------------------------------------!
           !     local-orbital-local-orbital contribution   !
           !------------------------------------------------!
	   do j=1, 3
	      do ilo1=1, nlorb(is)
		 l1=lorbl(ilo1, is)
		 do m1=-l1, l1
		    lm1=idxlm(l1, m1)
		    zv2(:)=zzero
		    do ilo2=1, nlorb(is)
		       l3=lorbl(ilo2, is)
		       do m3=-l3, l3
			  lm3=idxlm(l3, m3)
			  call zaxpy(nstfv, &
			       zone * riplolo(ilo1, m1, ilo2, m3, ias, j), &
			       locmt(1, ilo2, m3, ias), 1, zv2, 1)
		       end do
		    end do
		    call zoutpr(nstfv, nstfv, &
			 zone, locmt(1, ilo1, m1, ias), zv2, pm(1, 1, j))
		 end do
	      end do
	   end do
           ! end case of local orbitals
	end if
        ! end loop over atoms and species
     end do
  end do
  ! multiply y-component with imaginary unit
  pm(:, :, 2)=zi*pm(:, :, 2)
  !  calculate momentum matrix elements in the interstitial region
  forall (ist1=1:nstfv)
     evecfv1(ist1, :)=conjg(evecfv(1:ngp, ist1))
  end forall
  evecfv2(:, :)=evecfv(1:ngp, :)
  do j=1, 3
     do igp1=1, ngp
	ig1=igpig(igp1)
	iv1(:)=ivg(:, ig1)
	do igp2=1, ngp
	   ig2=igpig(igp2)
	   iv(:)=iv1(:)-ivg(:, ig2)
	   ig=ivgig(iv(1), iv(2), iv(3))
	   cfunt(igp1, igp2)=zi*vgpc(j, igp2)*cfunig(ig)
	end do
     end do
     call zgemm('n', 'n', ngp, nstfv, ngp, zone, cfunt, &
	  ngp, evecfv2, ngp, zzero, h, ngp)
     call zgemm('n', 'n', nstfv, nstfv, ngp, zone, evecfv1, &
	  nstfv, h, ngp, zzero, pmt, nstfv)
     pm(:, :, j)=pm(:, :, j)+pmt(:, :)
  end do
  ! multiply by -i and set lower triangular part
  do ist=1, nstfv
     do jst=ist, nstfv
	pm(ist, jst, :)=-zi*pm(ist, jst, :)
	pm(jst, ist, :)=conjg(pm(ist, jst, :))
     end do
  end do
  ! compute the second-variational momentum matrix elements
  if (input%groundstate%tevecsv) then
     do i=1, nstsv
	do j=1, nstsv
	   zv(:)=0.d0
	   k=0
	   do ispn=1, nspinor
	      do ist=1, nstfv
		 k=k+1
		 l=(ispn-1)*nstfv
		 do jst=1, nstfv
		    l=l+1
		    zt1=conjg(evecsv(k, i))*evecsv(l, j)
		    zv(:)=zv(:)+zt1*pm(ist, jst, :)
		 end do
	      end do
	   end do
	   pmat(:, i, j)=zv(:)
	end do
     end do
  else
     do j=1, 3
	pmat(j, :, :)=pm(:, :, j)
     end do
  end if
  deallocate(wfmt, gwfmt, pm, cfunt, h, pmt, evecfv1, evecfv2)
end subroutine genpmat2
!EOC
