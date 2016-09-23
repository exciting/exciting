!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: fxc_alda_check
! !INTERFACE:
!
!
Subroutine fxc_alda_check
! !USES:
      Use modinput
      Use modmain, Only: lmmaxvr, ngrtot, rhoir, rhomt, vxcir, vxcmt,xctype
      Use modxcifc, Only: xcifc
      Use modfxcifc
! !DESCRIPTION:
!   Checks the validity of the analytical expressions for the ALDA
!   exchange-correlation kernel.
!
! !REVISION HISTORY:
!   Created April 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Real (8), Allocatable :: vx (:), ex (:)
      Real (8), Allocatable :: vc (:), ec (:)
      Real (8), Allocatable :: dvx (:), dvc (:), vxc (:), dvxc2 (:), &
     & dvxc (:), cf (:, :)
      Integer :: m, nrho, irh
      Real (8), Allocatable :: rhogr (:)
      Real (8) :: rhoint (2)
!
      Call init0
      Call init1
      Call readstate
      Call init2
!
      m = Max (lmmaxvr, ngrtot)
!
      nrho = 100000
      rhoint (1) = 1.d0
      rhoint (2) = 7400.0d0
!
      Allocate (rhogr(nrho))
      Allocate (dvx(nrho), dvc(nrho), dvxc2(nrho), dvxc(nrho), cf(3, &
     & nrho))
      Allocate (ex(nrho), ec(nrho), vx(nrho), vc(nrho), vxc(nrho))
!
      Do irh = 1, nrho
         rhogr (irh) = rhoint (1) + (rhoint(2)-rhoint(1)) * dble &
        & (irh-1) / dble (nrho)
      End Do
!
      Call xcifc (xctype, n=nrho, rho=rhogr, &
     & ex=ex, ec=ec, vx=vx, vc=vc)
      vxc (:) = vx (:) + vc (:)
!
  ! use analytic expression of xc-kernel
      Call xcd_pwca (nrho, rhogr, dvx, dvc)
      dvxc (:) = dvx (:) + dvc (:)
!
  ! numerical differentiation
      Call fderiv (1, nrho, rhogr, vxc, dvxc2, cf)
!
  ! plot both versions
!
      Do irh = 1, nrho
         Write (2000, '(i6, 100g18.10)') irh, rhogr (irh), vxc (irh), &
        & dvxc2 (irh), dvxc (irh)
      End Do
!
      Write (*,*) 'minimum interstital density:', minval (rhoir)
      Write (*,*) 'minimum muffin-tin density :', minval (rhomt)
      Write (*,*) 'maximum interstital density:', maxval (rhoir)
      Write (*,*) 'maximum muffin-tin density :', maxval (rhomt)
      Write (*,*)
      Write (*,*) 'minimum interstital potential:', minval (vxcir)
      Write (*,*) 'minimum muffin-tin potential :', minval (vxcmt)
      Write (*,*) 'maximum interstital potential:', maxval (vxcir)
      Write (*,*) 'maximum muffin-tin potential :', maxval (vxcmt)
!
End Subroutine fxc_alda_check
!EOC
