
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: reformatdynamicalmatrices
! !INTERFACE:
subroutine reformatdynamicalmatrices
! !USES:
use modinput
Use mod_atoms
use mod_qpoint
use mod_constants
! !DESCRIPTION:
!   Collecting pieces of the dynamical matrices and assembling them in a nicer
!   way, such that $3\times 3$ matrices are displayed for each combination of
!   atoms and each $\bf{q}$ point.
!
! !REVISION HISTORY:
!   Created February 2010 (Sagmeister)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: n
! allocatable arrays
      Complex (8), Allocatable :: dyn (:, :)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: dynmatq(:,:,:,:,:,:,:)
! initialise universal variables
      Call init0
      Call init2
      n = 3 * natmtot
      allocate(dyn(3,n))
      Allocate (dynq(n, n, nqpt))
      Allocate (dynmatq(3,3,natmmax,nspecies,natmmax,nspecies,nqpt))
! read in the dynamical matrices
      Call readdyn (.false.,dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices
      call writedynamicalmatrices('DYNMAT.OUT',dynmatq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices including sumrule correction
      call writedynamicalmatrices('DYNMAT_SUMRULE.OUT',dynmatq)
! read in the dynamical matrices including symmetrization
      Call readdyn (.true.,dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices
      call writedynamicalmatrices('DYNMAT_SYM.OUT',dynmatq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices including sumrule correction
      call writedynamicalmatrices('DYNMAT_SYM_SUMRULE.OUT',dynmatq)
      Deallocate (dyn, dynq, dynmatq)
      Write (*, '("Info(reformatdynamicalmatrices): reformatted dynamical matrices")')
      Write (*, '(" including symmetrization and/or application of accoustic sumrule")')
      Write (*, '(" written to DYNMAT.OUT, DYNMAT_SYM.OUT, DYNMAT_SUMRULE.OUT and DYNMAT_SYM_SUMRULE.OUT")')
      Write (*,*)
End Subroutine
!EOC


subroutine reorderdynmat(dynq,dynmatq)
  Use mod_atoms
  use mod_qpoint
  implicit none
! arguments
  Complex (8), intent(in) :: dynq (3*natmtot,3*natmtot,nqpt)
  Complex (8), intent(out) :: dynmatq(3,3,natmmax,nspecies,natmmax,nspecies,nqpt)
! local variables
  Integer :: n, iq, i, j, is, ia, ip, js, ja, jp
  n = 3 * natmtot
! store dynamical matrices in a different way
  Do iq = 1, nqpt
     i = 0
     Do is = 1, nspecies
        Do ia = 1, natoms (is)
           Do ip = 1, 3
              i = i + 1
              j = 0
              Do js = 1, nspecies
                 Do ja = 1, natoms (js)
                    Do jp = 1, 3
                       j = j + 1
                       dynmatq(ip,jp,ia,is,ja,js,iq)=dynq (i, j, iq)
                    End Do
                 End Do
              End Do
           End Do
! end loops over atoms and species
        End Do
     End Do
! end loop over q-vectors
  End Do
end subroutine


!BOP
! !ROUTINE: writedynamicalmatrices
! !INTERFACE:
subroutine writedynamicalmatrices(fname,dynmq)
! !USES:
use modinput
Use mod_atoms
use mod_qpoint
use mod_constants
! !DESCRIPTION:
!   Write out the dynamical matrices to file {\tt DYNMAT.OUT}.
!
! !REVISION HISTORY:
!   Created February 2010 (Sagmeister)
!EOP
!BOC
      Implicit None
! arguments
      character(*), intent(in) :: fname
      complex(8), intent(in) :: dynmq(3,3,natmmax,nspecies,natmmax,nspecies,nqpt)
! local variables
      integer :: ip,ia,is,ias,ja,js,jas,iq
      ! write out dynamical matrix
      open(50,file=trim(fname),form="formatted",action="write",status="replace")
      Do iq = 1, nqpt
         write(50,'("q-point, vql: ",i8,3g18.10)') iq, vql(:,iq)
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
            ias=idxas(ia,is)
               Do js = 1, nspecies
                  Do ja = 1, natoms (js)
                     jas=idxas(ja,js)
      Write (50, '("First atom: is = ", I4, " (", A, "), ia = ", I4, &
       &"; second atom: js = ", I4, " (", A, "), ja = ", I4)') &
       & is, trim(input%structure%speciesarray(is)%species%chemicalSymbol),ia, &
       & js, trim(input%structure%speciesarray(js)%species%chemicalSymbol),ja
      write(50,'("First and second atom (overall): ias, jas = ",i6,", ", &
       &i6,"; dynamical matrix below")') ias, jas
                     do ip=1,3
                        write(50,'(2g18.10,"   ",2g18.10,"   ",2g18.10)') &
                          dynmq(ip,:,ia,is,ja,js,iq)
                     end do
                     write(50,*)
! end loops over atoms and species
                  End Do
               End Do
! end loops over atoms and species
            End Do
         End Do
         write(50,*)
! end loop over q points
      end do
      close(50)
End Subroutine
!EOC

