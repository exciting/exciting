
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

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
