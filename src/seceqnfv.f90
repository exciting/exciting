!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
!
!
Subroutine seceqnfv (nmatp, ngp, igpig, vgpc, apwalm, evalfv, evecfv)
  ! !USES:
      Use modinput
      Use modmain
      Use modfvsystem
!
  ! !INPUT/OUTPUT PARAMETERS:
  !   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
  !   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
  !   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
  !   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  !   Solves the secular equation,
  !   $$ (H-\epsilon O)b=0, $$
  !   for the all the first-variational states of the input $k$-point.
  !
  ! !REVISION HISTORY:
  !   Created March 2004 (JKD)
  !EOP
  !BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: nmatp
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Out) :: evalfv (nstfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
  ! local variables
      Type (evsystem) :: system
      Logical :: packed
      integer ik,jk,ig,iv(3),ist
      Complex (8), allocatable :: kinetic(:,:),vectors(:,:)

  ! allocatable arrays

!
  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!
!

      packed = input%groundstate%solver%packedmatrixstorage

!      if (input%groundstate%ValenceRelativity.eq.'lkh') then
       Call newsystem (system, packed, nmatp,(input%groundstate%ValenceRelativity.eq.'lkh'))
!      else
!       Call newsystem (system, packed, nmatp,.false.)
!      endif
      Call hamiltonandoverlapsetup (system, ngp, apwalm, igpig, vgpc)
!
  !------------------------------------!
  !     solve the secular equation     !
  !------------------------------------!
if (.false.) then
      do ik=1,system%h1%rank
        do jk=1,system%h1%rank
          !write(*,*) real(system%h1%za(jk,ik)),real(system%hamilton%za(jk,ik)),real(system%overlap%za(jk,ik))
          write(*,*) real(system%h1%za(jk,ik)),imag(system%h1%za(jk,ik)),abs(system%h1%za(jk,ik))
        enddo
        read(*,*)
      enddo
endif
      if (input%groundstate%ValenceRelativity.eq.'lkh') then
        if (packed) then
          system%overlap%zap=system%overlap%zap+system%h1%zap
          system%hamilton%zap=system%hamilton%zap+system%h1%zap*input%groundstate%energyref
        else
          system%overlap%za=system%overlap%za+system%h1%za
          system%hamilton%za=system%hamilton%za+system%h1%za*input%groundstate%energyref
        endif
      endif
      
      Call solvewithlapack(system,nstfv,evecfv,evalfv)

End Subroutine seceqnfv
!EOC
