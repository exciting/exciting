!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: dumpparams_mpiiter
! !INTERFACE:
!
!
Subroutine dumpparams_mpiiter (string, comment)
! !USES:
      Use modinput
      Use modmain
      Use sclcontroll
! !DESCRIPTION:
!   Writes out all input parameters which can be specified in the input file
!   {\tt exciting.in}.
!   Only show those array elements that are within a corresponding cutoff.
!   Trailling whitespaces in string expressions are trimmed.
!   This routine refers only to the parameters related to the MPI
!   parallelization and the iterative solver.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Character (*), Intent (In) :: string, comment
      Open (Unit=77, File=trim(string), Action='write', Position='appen&
     &d')
      Write (77,*)
      Write (77, '("! EXCITING version ", I1.1, ".", I1.1, ".", I3.3)') &
     & version
      Write (77, '(a)') trim (comment)
      Write (77,*)
      Write (77, '("tarpack")')
      Write (77,*) tarpack
      Write (77,*)
      Write (77, '("packedmatrixstorage")')
      Write (77,*) input%groundstate%solver%packedmatrixstorage
      Write (77,*)
      Write (77, '("epsarpack")')
      Write (77,*) input%groundstate%solver%epsarpack
      Write (77,*)
      Write (77, '("epsresid")')
      Write (77,*) epsresid
      Write (77,*)
      Write (77, '("maxncv")')
      Write (77,*) maxncv
      Write (77,*)
      Write (77, '("lowesteval")')
      Write (77,*) lowesteval
      Close (77)
End Subroutine dumpparams_mpiiter
!EOC
