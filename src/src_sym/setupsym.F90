


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine setupsym
  use modmain
use modinput
  use modsym
  implicit none
  if (allocated(sgmut)) deallocate(sgmut)
  allocate(sgmut(nsymcrys, nsymcrys))
  ! generate multiplication table
  call gensymcmut(input%structure%epslat, maxsymcrys, nsymcrys, symlat, lsplsymc, vtlsymc, sgmut, &
     abelsg)
  if (allocated(negenr)) deallocate(negenr)
  allocate(negenr(nsymcrys))
  if (allocated(genr)) deallocate(genr)
  allocate(genr(nsymcrys))
  if (allocated(orbgenr)) deallocate(orbgenr)
  allocate(orbgenr(nsymcrys, nsymcrys))
  ! find generators   
  call findsymgenr(nsymcrys, sgmut, ngenr, negenr, genr, orbgenr)
end subroutine setupsym
