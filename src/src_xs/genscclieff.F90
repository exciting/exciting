



! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine genscclieff(iqr, nmax, n, scieff)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: iqr, n, nmax
  complex(8), intent(out) :: scieff(nmax, nmax)
  ! local variables
  logical :: tq0
  complex(8), allocatable :: scrn(:, :), scrnw(:, :, :), scrnh(:, :)
  logical, external :: tqgamma
  allocate(scrn(n, n), scrnw(n, 2, 3), scrnh(3, 3))
  ! read screening from file
  call getscreen(iqr, n, scrnh, scrnw, scrn)
  tq0=tqgamma(iqr)
  if (tq0) then
     ! averaging using Lebedev-Laikov spherical grids
     call angavsc0(n, nmax, scrnh, scrnw, scrn, scieff)
  else
     ! averaging using numerical method and extrapolation
     call avscq(iqr, n, nmax, scrn, scieff)
  end if
end subroutine genscclieff


!//////////////////////////////////////////////////////////////////////////////


subroutine angavsc0(n, nmax, scrnh, scrnw, scrn, scieff)
  use modmain
use modinput
  use modxs
  use invert
  implicit none
  ! arguments
  integer, intent(in) :: n, nmax
  complex(8), intent(in) :: scrn(n, n), scrnw(n, 2, 3), scrnh(3, 3)
  complex(8), intent(out) :: scieff(nmax, nmax)
  ! local variables
  integer, parameter :: nsphcov=5810, iq0=1
  integer :: iop, jop, j1, j2, ji, itp, lm, ntpsph
  real(8) :: vomega, t00, r, qsz
  complex(8) :: dtns(3, 3), w1, w2
  real(8), allocatable :: plat(:, :), p(:), tp(:, :), spc(:, :), w(:)
  complex(8), allocatable :: m00lm(:), mx0lm(:), mxxlm(:)
  complex(8), allocatable :: ei00(:), eix0(:), ei0x(:), eixx(:)
  complex(8), allocatable :: ei00lm(:), eix0lm(:), ei0xlm(:), eixxlm(:)
  complex(8), allocatable :: ylm(:), zylm(:, :)
  complex(8), allocatable :: b(:, :), bi(:, :), u(:, :), v(:, :), s(:, :), t(:, :), e3(:, :), ie3(:, :)


  integer :: i1, i2
!!$  ! *** values for PA ***
!!$  call preset_dielten

  ! crystal volume
  vomega=omega*product(ngridq)
  ! Wigner-Seitz radius and spherical approximation to 1/q^2 average
  qsz=(6*pi**2/vomega)**(1.d0/3.d0)
  ! weight for sqrt(4pi)/q based on Wigner-Seitz radius
  w1=qsz**2*vomega/(4.d0*pi**2)*sqrt(fourpi)
  ! weight for 4pi/q^2 based on Wigner-Seitz radius
  w2=2*qsz*vomega/pi

  ! calculate RPA dielectric tensor including local field effects
  dielten0(:, :)=scrnh(:, :)
  if (n.gt.1) then
     allocate(b(n-1, n-1), bi(n-1, n-1), u(n-1, 3), v(3, n-1), s(n-1, 3), t(3, n-1))
     ! body of dielectric matrix
     b(:, :)=scrn(2:, 2:)
     ! column wing
     u(:, :)=scrnw(2:, 2, :)
     ! row wing
     v(:, :)=transpose(scrnw(2:, 1, :))
     select case(input%xs%BSE%scrherm)
     case(0)
        ! use full matrix (both wings and full body)
     case(1)
        ! Hermitian average matrix (average both wings and body)
	b=0.5d0*(b+conjg(transpose(b)))
	u=0.5d0*(u+conjg(transpose(v)))
	v=conjg(transpose(u))
     case(2)
        ! use upper triangle (take row wing, assign column wing)
	u=conjg(transpose(v))
     case(3)
        ! use lower triangle (take column wing, assign row wing)
	v=conjg(transpose(u))
     case default
	write(*, *)
	write(*, '("Error(angavsc0): not a valid flag:", i6)') input%xs%BSE%scrherm
	write(*, *)
	call terminate
     end select
     ! invert body (optionally including Hermitian average)
     call zinvert_hermitian(input%xs%BSE%scrherm, b, bi)
     s=matmul(bi, u)
     t=matmul(v, bi)
     dielten=dielten0-matmul(v, s)
  else
     dielten=dielten0
  end if
  ! symmetrize the dielectric tensor
  dtns(:, :)=dielten(:, :)
  do iop=1, 3
     do jop=1, 3
	call symt2app(iop, jop, 1, symt2, dtns, dielten(iop, jop))
     end do
  end do

  ! calculate averaged screened Coulomb interaction in Fourier space at Gamma point
  select case(trim(input%xs%BSE%sciavtype))
  case('spherical')
     ! scaling factor
     t00=(omega/(twopi)**3)*product(ngridq)
     ! number of points on sphere
     if (tleblaik) then
     ntpsph=input%xs%BSE%nleblaik
     else
	ntpsph=nsphcov
     end if
     if (lmmaxdielt.gt.ntpsph) then
	write(*, *)
	write(*, '("Error(angavdm0): lmmaxdielt.gt.ntpsph: ", 2i6)') lmmaxdielt, &
	     ntpsph
	write(*, *)
	stop
     end if
     allocate(plat(3, ntpsph), p(ntpsph))
     allocate(m00lm(lmmaxdielt), mx0lm(lmmaxdielt), mxxlm(lmmaxdielt))
     allocate(ei00(ntpsph), eix0(ntpsph), ei0x(ntpsph), eixx(ntpsph))
     allocate(ei00lm(lmmaxdielt), eix0lm(lmmaxdielt), ei0xlm(lmmaxdielt), eixxlm(lmmaxdielt))
     allocate(ylm(lmmaxdielt), zylm(ntpsph, lmmaxdielt))
     allocate(tp(2, ntpsph), spc(3, ntpsph))
     allocate(w(ntpsph))
     if (tleblaik) then
        ! generate Lebedev Laikov grid
	call leblaik(ntpsph, spc, w)
        ! generate tetha and phi angles
	do itp=1, ntpsph
	   call sphcrd(spc(:, itp), r, tp(:, itp))
	end do
     else
        ! distribution is assumed to be uniform
	w(:)=1.d0/ntpsph
        ! generate spherical covering set (angles and coordinates)
	call sphcover(ntpsph, tp)
	spc(1, :)=sin(tp(1, :))*cos(tp(2, :))
	spc(2, :)=sin(tp(1, :))*sin(tp(2, :))
	spc(3, :)=cos(tp(1, :))
     end if
     ! generate spherical harmonics on covering set
     do itp=1, ntpsph
     call genylm(input%xs%BSE%lmaxdielt, tp(:, itp), ylm)
	zylm(itp, :)=ylm(:)
     end do
     ! unit vectors of spherical covering set in lattice coordinates
     plat=matmul(binv, spc)
     ! distances to subcell cell boundaries in reciprocal space
     do itp=1, ntpsph
	p(itp:)=1.d0/(2.d0*maxval(abs(ngridq(:)*plat(:, itp)), 1))
     end do
     ! calculate function on covering set
     do itp=1, ntpsph
        ! head, 1/(p*L*p)
	ei00(itp)=1.d0/dot_product(spc(:, itp), matmul(dielten, spc(:, itp)))
     end do
     ! calculate lm-expansion coefficients
     do lm=1, lmmaxdielt
	ei00lm(lm)=fourpi*dot_product(zylm(:, lm), ei00*w)
	m00lm(lm)=fourpi*dot_product(zylm(:, lm), p*w)
	mx0lm(lm)=fourpi*dot_product(zylm(:, lm), p**2/2.d0*w)
	mxxlm(lm)=fourpi*dot_product(zylm(:, lm), p**3/3.d0*w)
     end do
     ! subcell average (head)
     scieff(1, 1)=fourpi*t00*dot_product(m00lm, ei00lm)
     ! loop over (G,Gp) indices
     do j1=2, n
	do itp=1, ntpsph
           ! wing, -p*S/(p*L*p)
	   eix0(itp)=-dot_product(spc(:, itp), s(j1-1, :))*ei00(itp)
           ! wing, -p*T/(p*L*p)
	   if (input%xs%BSE%scrherm.eq.0) ei0x(itp)=-dot_product(spc(:, itp), t(:, j1-1))*ei00(itp)
	end do
	do lm=1, lmmaxdielt
	   eix0lm(lm)=fourpi*dot_product(zylm(:, lm), eix0*w)
	   if (input%xs%BSE%scrherm.eq.0) ei0xlm(lm)=fourpi*dot_product(zylm(:, lm), ei0x*w)
	end do
        ! subcell average (wings)
	scieff(j1, 1)=sqrt(fourpi)*sptclg(j1, iq0)*t00*dot_product(mx0lm, eix0lm)
	if (input%xs%BSE%scrherm.eq.0) then
	  scieff(1, j1)=sqrt(fourpi)*sptclg(j1, iq0)*t00*dot_product(mx0lm, ei0xlm)
	else
	  scieff(1, j1)=conjg(scieff(j1, 1))
	end if
	if (input%xs%BSE%sciavbd) then
	   ji=j1
	   if (input%xs%BSE%scrherm.eq.0) ji=2
	   do j2=ji, n
	      do itp=1, ntpsph
                 ! body, B^-1 + p*S p*T/(p*L*p)
		 eixx(itp) = bi(j1 - 1, j2 - 1) + dot_product(spc(:, itp), s(j1 - 1, :))* &
		      dot_product(spc(:, itp), t(:, j2 - 1)) * ei00(itp)
	      end do
	      do lm=1, lmmaxdielt
		 eixxlm(lm)=fourpi*dot_product(zylm(:, lm), eixx*w)
	      end do
              ! subcell average (body)
	      scieff(j1, j2) = sptclg(j1, iq0) * sptclg(j2, iq0) * t00* &
		   dot_product(mxxlm, eixxlm)
	      if (input%xs%BSE%scrherm.ne.0) scieff(j2, j1)=conjg(scieff(j1, j2))
	   end do
	else
           ! no subcell average (body)
	   scieff(j1, 2:n)=bi(j1-1, :)
	end if
     end do
     deallocate(ei00, eix0, ei0x, eixx, ei00lm, eix0lm, ei0xlm, m00lm, mx0lm, mxxlm)
     deallocate(ylm, zylm, tp, spc, w, plat, p)
  case('screendiag', 'invscreendiag')
     if (input%xs%BSE%sciavbd) then
	write(*, *)
	write(*, '("Error(angavsc0): (inv)screendiag-method does not allow for averaging the body of W")')
	write(*, *)
	stop
     end if
     allocate(e3(n+2, n+2), ie3(n+2, n+2))
     ! invert dielectric matrix including 3 times G=0 according to the limits
     ! q->0_x, q->0_y, q->0_z
     ! G=0, G'=0 elements
     e3(1:3, 1:3)=dielten0(:, :)
     ! G!=0, G'=0 components and vice versa
     if (n.gt.1) then
	do i1=1, 3
	   do j2=2, n
	      e3(i1, j2+2)=scrnw(j2, 1, i1)
	   end do
	end do
	do j1=2, n
	   do i2=1, 3
	      e3(j1+2, i2)=scrnw(j1, 2, i2)
	   end do
	end do
	do j1=2, n
	   do j2=2, n
	      e3(j1+2, j2+2)=scrn(j1, j2)
	   end do
	end do
     end if
     call zinvert_hermitian(input%xs%BSE%scrherm, e3, ie3)
     ! select again
     select case(trim(input%xs%BSE%sciavtype))
     case('screendiag')
        ! head
	scieff(1, 1)=w2*1.d0/((e3(1, 1)+ie3(2, 2)+ie3(3, 3))/3.d0)
	if (n.gt.1) then
           ! wings, set to zero in this approximation
	   scieff(1, 2:n)=zzero
	   scieff(2:n, 1)=zzero
           ! body, only diagonal is assigned
	   scieff(2:n, 2:n)=zzero
	   forall (j1=2:n)
	      scieff(j1, j1)=sptclg(j1, iq0)**2/e3(j1+2, j1+2)
	   end forall
	end if
     case('invscreendiag')
        ! head
	scieff(1, 1)=w2*(ie3(1, 1)+ie3(2, 2)+ie3(3, 3))/3.d0
        ! wings
	if (n.gt.1) then
	   forall (j1=2:n)
	      scieff(j1, 1)=w1*sptclg(j1, iq0)*(ie3(j1+2, 1)+ie3(j1+2, 2)+ie3(j1+3, 3))/3.d0
	      scieff(1, j1)=conjg(scieff(j1, 1))
	   end forall
           ! body
	   forall (j1=2:n, j2=2:n)
	      scieff(j1, j2)=sptclg(j1, iq0)*sptclg(j2, iq0)*ie3(j1+2, j2+2)
	   end forall
	end if
     end select
     deallocate(e3, ie3)
  case('none')
     iop=1 !!!only x-component here!!!
     ! longitudinal treatment, three components of vanishing q (direction)
     allocate(e3(n, n), ie3(n, n))
     e3(1, 1)=scrnh(iop, iop)
     e3(1, 2:n)=scrnw(2:, 1, iop)
     e3(2:n, 1)=scrnw(2:, 2, iop)
     e3(2:n, 2:n)=scrn(2:, 2:)
     call zinvert_hermitian(input%xs%BSE%scrherm, e3, ie3)
     write(*, *) 'eps^{-1}_{00}=', ie3(1, 1)
     ! head
     scieff(1, 1)=w2*ie3(1, 1)
     ! wings
     if (n.gt.1) then
	forall (j1=2:n)
	   scieff(j1, 1)=w1*sptclg(j1, iq0)*ie3(j1, 1)
	   scieff(1, j1)=conjg(scieff(j1, 1))
	end forall
     	! body
	forall (j1=2:n, j2=2:n)
	   scieff(j1, j2)=sptclg(j1, iq0)*sptclg(j2, iq0)*ie3(j1+2, j2+2)
	end forall
     end if
     deallocate(e3, ie3)
  case default
     write(*, *)
     write(*, '("Error(angavsc0): invalid averaging method")')
     write(*, *)
     stop
  end select

  if (n.gt.1) deallocate(b, bi, u, v, s, t)

  call writedielt('DIELTENS', 1, 0.d0, dielten, 1)
  call writedielt('DIELTENS_NOSYM', 1, 0.d0, dtns, 1)

end subroutine angavsc0

!//////////////////////////////////////////////////////////////////////////////


subroutine avscq(iqr, n, nmax, scrn, scieff)
  use modmain
use modinput
  use modxs
  use invert
  implicit none
  ! arguments
  integer, intent(in) :: iqr, n, nmax
  complex(8), intent(in) :: scrn(n, n)
  complex(8), intent(out) :: scieff(nmax, nmax)
  ! local variables
  integer :: iqrnr, j1, j2, flg
  real(8) :: clwt
  ! find reduced q-point in non-reduced set
  iqrnr=iqmap(ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))
  ! invert dielectric matrix
  call zinvert_hermitian(input%xs%BSE%scrherm, scrn, scieff(:n, :n))
  do j1=1, n
     do j2=1, j1
	if ((input%xs%BSE%sciavqhd.and.(j1.eq.1).and.(j2.eq.1)).or. &
	     (input%xs%BSE%sciavqwg.and.(j1.ne.1).and.(j2.eq.1)).or. &
	     (input%xs%BSE%sciavqwg.and.(j1.eq.1).and.(j2.ne.1)).or. &
	     (input%xs%BSE%sciavqbd.and.(j1.ne.1).and.(j2.ne.1))) then
           ! numerical averaging on grids with extrapolation to continuum
	   flg=2
	else
           ! analytic expression, no averaging
	   flg=0
	end if
        ! generate the (averaged) symmetrized Coulomb potential
	call genwiqggp(flg, iqrnr, j1, j2, clwt)
        ! multiply with averaged Coulomb potential
	scieff(j1, j2)=scieff(j1, j2)*clwt
        ! set upper triangle
	scieff(j2, j1)=conjg(scieff(j1, j2))
     end do
  end do

end subroutine avscq


!//////////////////////////////////////////////////////////////////////////////


subroutine preset_dielten
  use modmain
  use modxs
  implicit none
  ! TODO: testing
!!$  real(8) :: r(3,3)
  ! preset dielectric tensor for testing
  dielten(:, :)=zzero
!!$!  (values are for trans-polyacetylene) from 2x2x16 k-point grid
!!$  dielten(1,:)=(/ 2.91911039, 0.00000000, 3.49765354 /)
!!$  dielten(2,:)=(/ 0.00000000, 2.79383654, 0.00000000 /)
!!$  dielten(3,:)=(/ 3.49765354, 0.00000000, 102.25001110 /)
!!$  dielten(1,:)=dielten(1,:)+zi*(/ 0.00000000, 0.00000000, 0.00000579 /)
!!$  dielten(2,:)=dielten(2,:)+zi*(/ 0.00000000, 0.00000000, 0.00000000 /)
!!$  dielten(3,:)=dielten(3,:)+zi*(/ -0.00000579,0.00000000, 0.00000000 /)
!  (values are for trans-polyacetylene), 4x4x32 k-point grid
  dielten(1, :)=(/ 2.91911039, 0.00000000, 0.49765354 /)
  dielten(2, :)=(/ 0.00000000, 2.79383654, 0.00000000 /)
  dielten(3, :)=(/ 0.49765354, 0.00000000, 52.725001110 /)
  dielten(1, :)=dielten(1, :)+zi*(/ 0.00000000, 0.00000000, 0.00000579 /)
  dielten(2, :)=dielten(2, :)+zi*(/ 0.00000000, 0.00000000, 0.00000000 /)
  dielten(3, :)=dielten(3, :)+zi*(/ -0.00000579, 0.00000000, 0.00000000 /)
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
