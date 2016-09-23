!
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Subroutine setupsym
      Use modmain
      Use modinput
      Use modsym
      Implicit None
      If (allocated(sgmut)) deallocate (sgmut)
      Allocate (sgmut(nsymcrys, nsymcrys))
  ! generate multiplication table
      Call gensymcmut (input%structure%epslat, maxsymcrys, nsymcrys, &
     & symlat, lsplsymc, vtlsymc, sgmut, abelsg, spainvsym)
      If (allocated(negenr)) deallocate (negenr)
      Allocate (negenr(nsymcrys))
      If (allocated(genr)) deallocate (genr)
      Allocate (genr(nsymcrys))
      If (allocated(orbgenr)) deallocate (orbgenr)
      Allocate (orbgenr(nsymcrys, nsymcrys))
  ! find generators
      Call findsymgenr (nsymcrys, sgmut, ngenr, negenr, genr, orbgenr)
End Subroutine setupsym
