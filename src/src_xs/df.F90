

! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: df
! !INTERFACE:


subroutine df
! !USES:
use modinput
  use modmain
  use modxs
  use modmpi
  use m_writegqpts
  use m_xsgauntgen
  use m_findgntn0
  use m_genfilname
! !DESCRIPTION:
!   Control routine for setting up the Kohn-Sham response function or the
!   microscopic dielectric function/matrix for all specified ${\bf q}$-points.
!   Can be run with MPI parallelization for ${\bf q}$-points.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  character(*), parameter :: thisnam='df'
  character(256) :: filex
  integer :: iq
  if (.not.tscreen) call genfilname(setfilext=.true.)
  call init0
  ! initialise universal variables
  call init1
  ! save Gamma-point variables
  call xssave0
  ! initialize q-point set
  call init2
  if (tscreen) then
     ! generate Gaunt coefficients
     call xsgauntgen(max(input%groundstate%lmaxapw, lolmax), input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
     ! find indices for non-zero Gaunt coefficients
     call findgntn0(max(input%xs%lmaxapwwf, lolmax), max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
  end if
  ! read Fermi energy
  call readfermi
  ! w-point parallelization for dielectric function
  if (tscreen) then
     nwdf=1
     call genparidxran('q', nqpt)
  else
     call genparidxran('w', nwdf)
  end if
  ! set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
  input%xs%emattype=1
  ! write out q-points
  call writeqpts
  ! loop over q-points
  do iq=qpari, qparf
     call genfilname(iq=iq, fileext=filex)
     ! call for q-point
     if (.not.input%xs%gather) call dfq(iq)
     if (tscreen) call writegqpts(iq, filex)
     write(unitout, '(a, i8)') 'Info('//thisnam//'): Kohn Sahm response &
	  &function finished for q - point:', iq
  end do
  ! synchronize
  if (.not.input%xs%gather) call barrier
  if ((procs.gt.1).and.(rank.eq.0).and.(.not.tscreen)) call dfgather
  if (.not.input%xs%gather) call barrier
  write(unitout, '(a)') "Info("//trim(thisnam)//"): Kohn-Sham response &
       &function finished"
  if (input%xs%gather) then
     write(unitout, '(a)') "Info("//trim(thisnam)//"): gather option: &
	  &exiting program"
     call xsfinit
     call terminate
  end if
  if (.not.tscreen) call genfilname(setfilext=.true.)
  if (tscreen) call findgntn0_clear
end subroutine df
!EOC
