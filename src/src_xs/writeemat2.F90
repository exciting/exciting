
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeemat2
  use modmain
  use modxs
  use modmpi
  use m_ematq2
  use m_tdgauntgen
  use m_findgntn0
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'writeemat2'
  integer :: iq
  ! initialise universal variables
  call init0
  call init1
  ! initialize q-point set
  call init2xs
  ! generate index ranges for parallel execution
  if ((task.ge.300).or.(task.le.399)) call genparidxran('k')
  if ((task.ge.400).or.(task.le.499)) call genparidxran('q')
  ! write q-point set
  if (rank.eq.0) call writeqpts
  ! read Fermi energy from file
  call readfermi
  ! generate Gaunt coefficients
  call tdgauntgen(lmaxmax,lmaxemat,lmaxmax)
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(lmaxapwtd,lmaxapwtd,lmaxemat,tdgnt)
  write(unitout,'(a,3i8)') 'Info('//thisnam//'): Gaunt coefficients generated &
       &within lmax values:', lmaxmax,lmaxemat,lmaxmax
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of q-points: ',nqpt
  call flushifc(unitout)
  ! loop over q-points
  do iq=qpari,qparf
     ! call for q-point
     call ematq2(iq)
     if (iq-qpari+1.le.nqpt/procs) call barrier
     write(unitout,'(a,i8)') 'Info('//thisnam//'): matrix elements of the &
          &plane wave finished for q-point:',iq
     call flushifc(unitout)
  end do
  call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): matrix elements of the &
       &plane wave finished"
  call flushifc(unitout)
  call findgntn0_clear
end subroutine writeemat2
