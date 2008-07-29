
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine df
  use modmain
  use modxs
  use modmpi
  use m_xsgauntgen
  use m_findgntn0
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='df'
  integer :: iq
  if (.not.tscreen) call genfilname(setfilext=.true.)
  if (calledxs.eq.1) call init0
  ! initialise universal variables
  call init1
  ! save Gamma-point variables
  call xssave0
  ! initialize q-point set
  call init2xs
  if (tscreen) then
     ! generate Gaunt coefficients
     call xsgauntgen(max(lmaxapw,lolmax),lmaxemat,max(lmaxapw,lolmax))
     ! find indices for non-zero Gaunt coefficients
     call findgntn0(max(lmaxapwwf,lolmax),max(lmaxapwwf,lolmax),lmaxemat,xsgnt)
  end if
  ! read Fermi energy
  call readfermi
  ! w-point parallelization for dielectric function
  call genparidxran('w',nwdf)
  if (tscreen) then
     nwdf=1
     call genparidxran('q',nqpt)
  end if
  ! set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
  emattype=1
  ! loop over q-points
  do iq=qpari,qparf
     ! call for q-point
     if (.not.gather) call dfq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn Sahm response &
          &function finished for q-point:',iq
  end do
  ! synchronize
  if (.not.gather) call barrier
  if ((procs.gt.1).and.(rank.eq.0).and.(.not.tscreen)) call dfgather
  if (.not.gather) call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Kohn-Sham response &
       &function finished"
  if (gather) then
     write(unitout,'(a)') "Info("//trim(thisnam)//"): gather option: &
          &exiting program"
     call xsfinit
     call terminate
  end if
  if (.not.tscreen) call genfilname(setfilext=.true.)
  if (tscreen) call findgntn0_clear
end subroutine df
