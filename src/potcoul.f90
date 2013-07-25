!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: potcoul
! !INTERFACE:
!
!
Subroutine potcoul
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are converted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, lmax, l, m, lm
      Complex (8) zrho0
! allocatable arrays
      Real (8), Allocatable :: jlgr (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :)
      Complex (8), Allocatable :: zvclir (:)
      Complex (8), Allocatable :: zrho2mt(:, :, :),charge(:),potential(:),cf(:,:) !,zvclmt (:, :, :)
      real(8), allocatable :: pot1(:),pot2(:),ch1(:),ch2(:)
!      write(*,*) 'ttt'
      Allocate &
     & (jlgr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, &
     & ngvec, nspecies))
      Allocate (zrhomt(lmmaxvr, nrmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrmtmax, natmtot))
      Allocate (zvclir(ngrtot))
      Allocate (zrho2mt(lmmaxvr, nrmtmax, natmtot))
! 
if (.false.) then
      Allocate (charge(nrmtmax))
      Allocate (potential(nrmtmax))
      Allocate (pot1(nrmtmax))
      Allocate (pot2(nrmtmax))
      Allocate (ch1(nrmtmax))
      Allocate (ch2(nrmtmax))
endif      
      
      
      

      Allocate(cf(3,nrmtmax))
! convert real muffin-tin charge density to complex spherical harmonic expansion
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               Call rtozflm (input%groundstate%lmaxvr, rhomt(:, ir, &
              & ias), zrhomt(:, ir, ias))
            End Do
         End Do
      End Do
! store real interstitial charge density in complex array
      zrhoir (:) = rhoir (:)
! compute the required spherical Bessel functions
      lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
      Call genjlgpr (lmax, gc, jlgr)
! solve the complex Poisson's equation

      Call zpotcoul (nrmt, nrmtmax, spnrmax, spr, 1, gc, jlgr, ylmg, &
     & sfacg, spzn, zrhomt, zrhoir, zvclmt, zvclir, zrho0)
!     stop
!Andris
!generate the density from the potential
if (.false.) then
      write(*,*) 'howdy'
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            write(*,*) 'atom', ias
            lm=0
            do l=0,input%groundstate%lmaxvr
              do m=-l,l
                lm=lm+1
                potential(1:nrmt(is))=zvclmt(lm,1:nrmt(is),ias)
                do ir=1,nrmt(is)
                   potential(ir)=spr(ir,is)*potential(ir)
                enddo
                pot1(:)=dble(potential(:))
                pot2(:)=dble(potential(:))

                Call fderiv (2, nrmt(is), spr(:,is), pot1, ch1 , cf)                
                Call fderiv (2, nrmt(is), spr(:,is), pot2, ch2 , cf)
                charge(:)=dcmplx(ch1(:),ch2(:))
                do ir=1,nrmt(is)
                   charge(ir)=-(charge(ir)-dble(l*(l+1))*potential(ir)/spr(ir,is)**2)/spr(ir,is)*y00**2
                enddo
                write(*,*) lm, ias
                do ir=1,nrmt(is)
                 write(*,*) dble(charge(ir)),dble(zrhomt(lm, ir, ias))
!                 write(*,*) spr(ir,is),dble(potential(ir)/spr(ir,is))
                enddo
!                read(*,*) 
              enddo
            enddo
         End Do
      End Do
!endif
      deallocate(charge,potential,pot1,pot2,ch1,ch2,cf)
endif
!      stop
! convert complex muffin-tin potential to real spherical harmonic expansion
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
!            write(*,*) zvclmt(:,nrmt(is),ias)
            Do ir = 1, nrmt (is)
               Call ztorflm (input%groundstate%lmaxvr, zvclmt(:, ir, &
              & ias), vclmt(:, ir, ias))
            End Do
         End Do
      End Do
! store complex interstitial potential in real array
      vclir (:) = dble (zvclir(:))
      Deallocate (jlgr, zrhomt, zrhoir, zvclmt, zvclir)
!      write(*,*) 'coulomb done'
      Return
End Subroutine
!EOC
