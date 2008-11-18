
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscclieff(iqr,nmax,n,scieff)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: iqr,n,nmax
  complex(8), intent(out) :: scieff(nmax,nmax)
  ! local variables
  logical :: tq0
  complex(8), allocatable :: scrn(:,:),scrnw(:,:,:),scrnh(:,:)
  logical, external :: tqgamma
  allocate(scrn(n,n),scrnw(n,2,3),scrnh(3,3))
  ! read screening from file
  call getscreen(iqr,n,scrnh,scrnw,scrn)
  tq0=tqgamma(iqr)
  if (tq0) then
     ! averaging using Lebedev-Laikov spherical grids
     call angavsc0(n,nmax,scrnh,scrnw,scrn,scieff)
  else
     ! averaging using numerical method and extrapolation
     call avscq(iqr,n,nmax,scrn,scieff)
  end if
end subroutine genscclieff


!//////////////////////////////////////////////////////////////////////////////


subroutine avscq(iqr,n,nmax,scrn,scieff)
  use modmain
  use modxs
  use invert
  implicit none
  ! arguments
  integer, intent(in) :: iqr,n,nmax
  complex(8), intent(in) :: scrn(n,n)
  complex(8), intent(out) :: scieff(nmax,nmax)
  ! local variables
  integer :: iqrnr,j1,j2,flg
  real(8) :: clwt
  ! find reduced q-point in non-reduced set
  iqrnr=iqmap(ivqr(1,iqr),ivqr(2,iqr),ivqr(3,iqr))
  ! invert dielectric matrix
  call zinvert_hermitian(scrherm,scrn,scieff(:n,:n))
  do j1=1,n
     do j2=1,j1
        if ((sciavqhd.and.(j1.eq.1).and.(j2.eq.1)).or. &
             (sciavqwg.and.(j1.ne.1).and.(j2.eq.1)).or. &
             (sciavqbd.and.(j1.ne.1).and.(j2.ne.1))) then
           ! numerical averaging on grids with extrapolation to continuum
           flg=2
        else
           ! analytic expression, no averaging
           flg=0
        end if
        ! generate the (averaged) symmetrized Coulomb potential
        call genwiqggp(flg,iqrnr,j1,j2,clwt)
        ! multiply with averaged Coulomb potential
        scieff(j1,j2)=scieff(j1,j2)*clwt
        ! set upper triangle
        scieff(j2,j1)=conjg(scieff(j1,j2))
	
	if (abs(scieff(j1,j2)).gt.1.d5) then
	  write(*,'(a,3i5,2g18.10)') &
	   'scieff,iqr,j1,j2',iqr,j1,j2,abs(scieff(j1,j2))
	end if
	
     end do
  end do
end subroutine avscq


!//////////////////////////////////////////////////////////////////////////////


subroutine angavsc0(n,nmax,scrnh,scrnw,scrn,scieff)
  use modmain
  use modxs
  use invert
  implicit none
  ! arguments
  integer, intent(in) :: n,nmax
  complex(8), intent(in) :: scrn(n,n),scrnw(n,2,3),scrnh(3,3)
  complex(8), intent(out) :: scieff(nmax,nmax)
  ! local variables
  integer, parameter :: nsphcov=5810,iq0=1
  integer :: j1,j2,itp,lm,ntpsph
  real(8) :: t00,r
  real(8), allocatable :: plat(:,:),p(:),tp(:,:),spc(:,:),w(:)
  complex(8), allocatable :: m00lm(:),mx0lm(:),mxxlm(:)
  complex(8), allocatable :: ei00(:),eix0(:),eixx(:)
  complex(8), allocatable :: ei00lm(:),eix0lm(:),eixxlm(:)
  complex(8), allocatable :: ylm(:),zylm(:,:)
  complex(8), allocatable :: b(:,:),bi(:,:),u(:,:),s(:,:)
  ! scaling factor
  t00=(omega/(twopi)**3)*product(ngridq)
  
!!$  ! *** values for PA ***
!!$  call preset_dielten
!!$  zt1=1.d0/((dielten(1,1)+dielten(2,2)+dielten(3,3))/3.d0),
!!$  zt2=(1.d0/dielten(1,1)+1.d0/dielten(2,2)+1.d0/dielten(3,3))/3.d0
!!$  ! Wigner-Seitz radius and spherical approximation to 1/q^2 average
!!$  qsz=(6*pi**2/(omega*product(ngridq)))**(1.d0/3.d0)
!!$  clwt=2*qsz*omega*product(ngridq)/pi

  ! invert dielectric tensor
  dielten0(:,:)=scrnh(:,:)
  if (n.gt.1) then
     allocate(b(n-1,n-1),bi(n-1,n-1),u(n-1,3),s(n-1,3))
     ! body of dielectric matrix
     b(:,:)=scrn(2:,2:)
     ! wings of dielectric matrix
     u(:,:)=conjg(scrnw(2:,1,:))
     ! invert body (optionally including Hermitian average)
     call zinvert_hermitian(scrherm,b,bi)
     s=matmul(bi,u)
     dielten=dielten0-matmul(conjg(transpose(u)),s)
  else
     dielten=dielten0
  end if
  ! number of points on sphere
  if (tleblaik) then
     ntpsph=nleblaik
  else
     ntpsph=nsphcov
  end if
  if (lmmaxdielt.gt.ntpsph) then
     write(*,*)
     write(*,'("Error(angavdm0): lmmaxdielt.gt.ntpsph: ",2i6)') lmmaxdielt, &
          ntpsph
     write(*,*)
     stop
  end if
  allocate(plat(3,ntpsph),p(ntpsph))
  allocate(m00lm(lmmaxdielt),mx0lm(lmmaxdielt),mxxlm(lmmaxdielt))
  allocate(ei00(ntpsph),eix0(ntpsph),eixx(ntpsph))
  allocate(ei00lm(lmmaxdielt),eix0lm(lmmaxdielt),eixxlm(lmmaxdielt))
  allocate(ylm(lmmaxdielt),zylm(ntpsph,lmmaxdielt))
  allocate(tp(2,ntpsph),spc(3,ntpsph))
  allocate(w(ntpsph))
  if (tleblaik) then
     ! generate Lebedev Laikov grid
     call leblaik(ntpsph,spc,w)
     ! generate tetha and phi angles
     do itp=1,ntpsph
        call sphcrd(spc(:,itp),r,tp(:,itp))
     end do
  else
     ! distribution is assumed to be uniform
     w(:)=1.d0/ntpsph
     ! generate spherical covering set (angles and coordinates)
     call sphcover(ntpsph,tp)
     spc(1,:)=sin(tp(1,:))*cos(tp(2,:))
     spc(2,:)=sin(tp(1,:))*sin(tp(2,:))
     spc(3,:)=cos(tp(1,:))
  end if
  ! generate spherical harmonics on covering set
  do itp=1,ntpsph
     call genylm(lmaxdielt,tp(:,itp),ylm)
     zylm(itp,:)=ylm(:)
  end do
  ! unit vectors of spherical covering set in lattice coordinates
  plat=matmul(binv,spc)
  ! distances to subcell cell boundaries in reciprocal space
  do itp=1,ntpsph
     p(itp:)=1.d0/(2.d0*maxval(abs(ngridq(:)*plat(:,itp)),1))
  end do
  ! calculate function on covering set
  do itp=1,ntpsph
     ! head, 1/(p*L*p)
     ei00(itp)=1.d0/dot_product(spc(:,itp),matmul(dielten,spc(:,itp)))
  end do
  ! calculate lm-expansion coefficients
  do lm=1,lmmaxdielt
     ei00lm(lm)=fourpi*dot_product(zylm(:,lm),ei00*w)
     m00lm(lm)=fourpi*dot_product(zylm(:,lm),p*w)
     mx0lm(lm)=fourpi*dot_product(zylm(:,lm),p**2/2.d0*w)
     mxxlm(lm)=fourpi*dot_product(zylm(:,lm),p**3/3.d0*w)
  end do
  ! subcell average (head)
  scieff(1,1)=fourpi*t00*dot_product(m00lm,ei00lm)
  ! loop over (G,Gp) indices
  do j1=2,n
     do itp=1,ntpsph
        ! wing, -p*S/(p*L*p)
        eix0(itp)=-dot_product(spc(:,itp),s(j1-1,:))*ei00(itp)
     end do
     do lm=1,lmmaxdielt
        eix0lm(lm)=fourpi*dot_product(zylm(:,lm),eix0*w)
     end do
     ! subcell average (wings)
     scieff(j1,1)=sqrt(fourpi)*sptclg(j1,iq0)*t00*dot_product(mx0lm,eix0lm)
     scieff(1,j1)=conjg(scieff(j1,1))
     if (sciavbd) then
        do j2=j1,n
           do itp=1,ntpsph
              ! body, B^-1 + p*S p*conjg(S)/(p*L*p)
              eixx(itp)=bi(j1-1,j2-1)*dot_product(spc(:,itp),s(j1-1,:))* &
                   dot_product(s(j2-1,:),spc(:,itp))*ei00(itp)
           end do
           do lm=1,lmmaxdielt
              eixxlm(lm)=fourpi*dot_product(zylm(:,lm),eixx*w)
           end do
           ! subcell average (body)
           scieff(j1,j2)=sptclg(j1,iq0)*sptclg(j2,iq0)*t00* &
                dot_product(mxxlm,eixxlm)
           scieff(j2,j1)=conjg(scieff(j1,j2))
        end do
     else
        ! no subcell average (body)
        scieff(j1,2:)=bi(j1-1,:)
     end if
  end do
  deallocate(ei00,ei00lm,m00lm,mx0lm,mxxlm,ylm,zylm,tp,spc,w,plat,p)
end subroutine angavsc0


!//////////////////////////////////////////////////////////////////////////////

subroutine preset_dielten
  use modmain
  use modxs
  implicit none
  real(8) :: r(3,3)
  ! preset dielectric tensor for testing
  dielten(:,:)=zzero
!!$!  (values are for trans-polyacetylene) from 2x2x16 k-point grid
!!$  dielten(1,:)=(/ 2.91911039, 0.00000000, 3.49765354 /)
!!$  dielten(2,:)=(/ 0.00000000, 2.79383654, 0.00000000 /)
!!$  dielten(3,:)=(/ 3.49765354, 0.00000000, 102.25001110 /)
!!$  dielten(1,:)=dielten(1,:)+zi*(/ 0.00000000, 0.00000000, 0.00000579 /)
!!$  dielten(2,:)=dielten(2,:)+zi*(/ 0.00000000, 0.00000000, 0.00000000 /)
!!$  dielten(3,:)=dielten(3,:)+zi*(/ -0.00000579,0.00000000, 0.00000000 /)
!  (values are for trans-polyacetylene), 4x4x32 k-point grid
  dielten(1,:)=(/ 2.91911039, 0.00000000, 0.49765354 /)
  dielten(2,:)=(/ 0.00000000, 2.79383654, 0.00000000 /)
  dielten(3,:)=(/ 0.49765354, 0.00000000, 52.725001110 /)
  dielten(1,:)=dielten(1,:)+zi*(/ 0.00000000, 0.00000000, 0.00000579 /)
  dielten(2,:)=dielten(2,:)+zi*(/ 0.00000000, 0.00000000, 0.00000000 /)
  dielten(3,:)=dielten(3,:)+zi*(/ -0.00000579,0.00000000, 0.00000000 /)
!!$! testing only
!!$  dielten(1,:)=(/ 3.0, 0.0, 0.0 /)
!!$  dielten(2,:)=(/ 0.0, 3.0, 0.0 /)
!!$  dielten(3,:)=(/ 0.0, 0.0, 3.0 /)
!!$  dielten(1,:)=dielten(1,:)+zi*(/ 0.0, 0.0, 0.0 /)
!!$  dielten(2,:)=dielten(2,:)+zi*(/ 0.0, 0.0, 0.0 /)
!!$  dielten(3,:)=dielten(3,:)+zi*(/ 0.0, 0.0, 0.0 /)
!!$  call random_number(r)
!!$  dielten(:,:)=dielten(:,:)+r(:,:)*1.d0
!!$  call random_number(r)
!!$  dielten(:,:)=dielten(:,:)+zi*r(:,:)*1.d0
end subroutine preset_dielten
