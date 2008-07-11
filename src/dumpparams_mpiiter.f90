
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dumpparams_mpiiter
! !INTERFACE:
subroutine dumpparams_mpiiter(string,comment)
! !USES:
  use modmain
  use sclcontroll
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
  implicit none
  ! arguments
  character(*), intent(in) :: string,comment
  open(unit=77,file=trim(string),action='write',status='replace')
  write(77,*)
  write(77,'("! EXCITING version ",I1.1,".",I1.1,".",I3.3)') version
  write(77,'(a)') trim(comment)
  write(77,*)
  write(77,'("iterativetype")')
  write(77,*) iterativetype
  write(77,*)
  write(77,'("packedmatrixstorage")')
  write(77,*) packedmatrixstorage
  write(77,*)
  write(77,'("epsarpack")')
  write(77,*) epsarpack
  write(77,*)
  write(77,'("epsresid")')
  write(77,*) epsresid
  write(77,*)
  write(77,'("maxncv")')
  write(77,*) maxncv
  write(77,*)
  write(77,'("lowesteval")')
  write(77,*) lowesteval
  close(77)
end subroutine dumpparams_mpiiter
!EOC
