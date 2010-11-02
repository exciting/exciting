

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_sp_2a
! !INTERFACE:


subroutine ggamt_sp_2a(is, rhoup, rhodn, g2up, g2dn, gvup, gvdn, gup2, gdn2, gupdn)
! !USES:
use modinput
use mod_Gvector
use mod_muffin_tin
use mod_SHT
use mod_atoms
use mod_potential_and_density
! !DESCRIPTION:
!   Computes the muffin-tin gradients $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho^{\uparrow}$,
!   $\nabla\rho^{\downarrow}$, $(\nabla\rho^{\uparrow})^2$,
!   $(\nabla\rho^{\downarrow})^2$ and
!   $\nabla\rho^{\uparrow}\cdot\nabla\rho^{\downarrow}$, which are passed in to
!   GGA functional subroutines of type 2. The exchange-correlation energy in
!   these routines has the functional form
!   $$ E_{xc}[\rho^{\uparrow},\rho^{\downarrow}]=\int d^3r\,\hat{\epsilon}_{xc}
!    \bigl(\rho^{\uparrow}({\bf r}),\rho^{\downarrow}({\bf r}),
!    (\nabla\rho^{\uparrow}({\bf r}))^2,(\nabla\rho^{\downarrow}({\bf r}))^2,
!    \nabla\rho^{\uparrow}({\bf r})
!    \cdot\nabla\rho^{\downarrow}({\bf r})\bigr), $$
!   where $\hat{\epsilon}_{xc}({\bf r})=\epsilon_{xc}({\bf r})\rho({\bf r})$ is
!   the xc energy per unit volume, with $\epsilon_{xc}$ being the xc energy per
!   electron, and $\rho=\rho^{\uparrow}+\rho^{\downarrow}$. From the gradients
!   above, type 2 GGA routines return $\epsilon_{xc}$, but not directly the xc
!   potentials. Instead they generate the derivatives
!   $\partial\hat{\epsilon}_{xc}/\partial\rho^{\uparrow}({\bf r})$,
!   $\partial\hat{\epsilon}_{xc}/\partial(\nabla\rho^{\uparrow}({\bf r}))^2$,
!   and the same for down spin, as well as
!   $\partial\hat{\epsilon}_{xc}/\partial(\nabla\rho^{\uparrow}({\bf r})
!   \cdot\nabla\rho^{\downarrow}({\bf r}))$. In a post-processing step invoked
!   by {\tt ggamt\_sp\_2b}, integration by parts is used to obtain the xc
!   potential explicitly with
!   \begin{align*}
!    V_{xc}^{\uparrow}({\bf r})=&\frac{\partial\hat{\epsilon}_{xc}}{\partial
!    \rho^{\uparrow}({\bf r})}-2\left(\nabla\frac{\partial\hat{\epsilon}_{xc}}
!    {\partial(\nabla\rho^{\uparrow})^2}\right)\cdot\nabla\rho^{\uparrow}
!    -2\frac{\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow})^2}\nabla^2
!    \rho^{\uparrow}\&
!    &-\left(\nabla\frac{\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow}
!    \cdot\nabla\rho^{\downarrow})}\right)\cdot\nabla\rho^{\downarrow}
!    -\frac{\partial\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow}\cdot
!    \nabla\rho^{\downarrow})}\nabla^2\rho^{\downarrow},
!   \end{align*}
!   and similarly for $V_{xc}^{\downarrow}$.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
real(8), intent(in) :: rhoup(lmmaxvr, nrmtmax)
real(8), intent(in) :: rhodn(lmmaxvr, nrmtmax)
real(8), intent(in) :: g2up(lmmaxvr, nrmtmax)
real(8), intent(in) :: g2dn(lmmaxvr, nrmtmax)
real(8), intent(out) :: gvup(lmmaxvr, nrmtmax, 3)
real(8), intent(out) :: gvdn(lmmaxvr, nrmtmax, 3)
real(8), intent(out) :: gup2(lmmaxvr, nrmtmax)
real(8), intent(out) :: gdn2(lmmaxvr, nrmtmax)
real(8), intent(out) :: gupdn(lmmaxvr, nrmtmax)
! local variables
integer::nr, i
! allocatable arrays
real(8), allocatable :: rfmt1(:, :), rfmt2(:, :), grfmt(:, :, :)
allocate(rfmt1(lmmaxvr, nrmtmax), rfmt2(lmmaxvr, nrmtmax))
allocate(grfmt(lmmaxvr, nrmtmax, 3))
nr=nrmt(is)
!----------------!
!     rho up     !
!----------------!
! convert rhoup to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, rhoup, lmmaxvr, 0.d0, &
 rfmt1, lmmaxvr)
! compute grad^2 rhoup in spherical coordinates
call grad2rfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, rfmt1, rfmt2)
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, rfmt2, lmmaxvr, 0.d0, &
 g2up, lmmaxvr)
! grad rhoup in spherical coordinates
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, gvup(:, :, i), lmmaxvr)
end do
! (grad rhoup)^2
gup2(:, 1:nr)=gvup(:, 1:nr, 1)**2+gvup(:, 1:nr, 2)**2+gvup(:, 1:nr, 3)**2
!------------------!
!     rho down     !
!------------------!
! convert rhodn to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, rhodn, lmmaxvr, 0.d0, &
 rfmt1, lmmaxvr)
! compute grad^2 rhodn in spherical coordinates
call grad2rfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, rfmt1, rfmt2)
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, rfmt2, lmmaxvr, 0.d0, &
 g2dn, lmmaxvr)
! grad rhodn in spherical coordinates
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, gvdn(:, :, i), lmmaxvr)
end do
! (grad rhodn)^2
gdn2(:, 1:nr)=gvdn(:, 1:nr, 1)**2+gvdn(:, 1:nr, 2)**2+gvdn(:, 1:nr, 3)**2
! (grad rhoup).(grad rhodn)
gupdn(:, 1:nr) = gvup(:, 1:nr, 1) * gvdn(:, 1:nr, 1) &
	     +gvup(:, 1:nr, 2) * gvdn(:, 1:nr, 2) &
	     +gvup(:, 1:nr, 3) * gvdn(:, 1:nr, 3)
deallocate(rfmt1, rfmt2, grfmt)
return
end subroutine
!EOC
