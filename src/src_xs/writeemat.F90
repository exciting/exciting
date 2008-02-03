
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeemat
  use modmain
  use modxs
  use modmpi
  use m_tdgauntgen
  use m_findgntn0
  use m_filedel
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='writeemat'
  integer :: iq
  ! initialise universal variables
  if (calledxs.eq.1) call init0
  call init1
  call init2xs
  ! k-point parallelization for TDDFT
  if ((task.ge.300).or.(task.le.399)) partype='k'
  ! q-point parallelization for screening
  if ((task.ge.400).or.(task.le.499)) partype='q'
  call genparidxran(partype)
   ! write q-point set
  if (rank.eq.0) call writeqpts
  ! read Fermi energy from file
  call readfermi
  ! save variables for the Gamma q-point
  call tdsave0
  ! generate Gaunt coefficients
  call tdgauntgen(lmaxapw,lmaxemat,lmaxapw)
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(lmaxapwtd,lmaxapwtd,lmaxemat,tdgnt)
  write(unitout,'(a,3i8)') 'Info('//thisnam//'): Gaunt coefficients generated &
       &within lmax values:', lmaxapw,lmaxemat,lmaxapw
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of q-points: ',nqpt
  call flushifc(unitout)
  if (gather) goto 10
  ! loop over q-points
  do iq=1,nqpt
     ! call for q-point
     call ematq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): matrix elements of the &
          &exponentials finished for q-point:',iq
     call flushifc(unitout)
  end do
  ! synchronize
  call barrier
10 continue
  ! gather from processes
  if ((procs.gt.1).and.(rank.eq.0).and.(partype.eq.'k')) then
     call ematgather
  end if
  ! synchronize
  call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): matrix elements of &
       &exponential expression finished"
  call findgntn0_clear
  call genfilname(setfilext=.true.)
end subroutine writeemat
