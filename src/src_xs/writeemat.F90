
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeemat
  use modmain
  use modtddft
  use modpar
  use m_ematq
  use m_tdgauntgen
  use m_findgntn0
  use m_getunit
  use m_filedel
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'writeemat'
  integer :: iq,ipar,un,qi
  logical :: existent

  ! initialise universal variables
  call init0
  call init1

  ! initialize q-point set
  call init2td()

  ! k-point interval for process
  call getrange(rank,nproc,nkpt,kpari,kparf)
  resumechkpts(1,2)=kpari
  resumechkpts(1,3)=kparf
  resumechkpts(2,2)=1
  resumechkpts(2,3)=nqpt

  ! write q-point set
  if (rank == 1) call writeqpts()

  ! read Fermi energy from file
  call readfermi

  ! save variables for the Gamma q-point
  call tdsave0()

  ! generate Gaunt coefficients
  call tdgauntgen(lmaxmax,lmaxemat,lmaxmax)

  ! find indices for non-zero Gaunt coefficients
  call findgntn0(lmaxapwtd,lmaxapwtd,lmaxemat,tdgnt)

  write(unitout,'(a,3i8)') 'Info('//thisnam//'): Gaunt coefficients generated &
       &within lmax values:', lmaxmax,lmaxemat,lmaxmax
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of q-points: ',nqpt

  if (gather) goto 10

  ! resume task, second checkpoint index is q-point index
  qi=1
  if (tresume) qi=resumechkpts(2,1)
  call getunit(un)
  ! loop over q-points
  do iq = qi, nqpt
     ! call for q-point
     call ematq(iq)
     resumechkpts(2,1)=iq
     call resupd(un,task,resumechkpts,' : q-point index')
     write(unitout,'(a,i8)') 'Info('//thisnam//'): matrix elements of the &
          &exponentials finished for q-point:',iq
  end do

  ! synchronize
  call getunit(un)
  call barrier(rank=rank,nproc=nproc,un=un,async=0,string='.barrier')

10 continue

  ! gather from processes
  if ((nproc.gt.1).and.(rank.eq.1)) call ematgather()
  if ((nproc.gt.1).and.(rank.eq.1)) call devalsvgather()

  ! synchronize
  call getunit(un)
  call barrier(rank=rank,nproc=nproc,un=un,async=0,string='.barrier')

  write(unitout,'(a)') "Info("//trim(thisnam)//"): matrix elements of &
       &exponential expression finished"

  call findgntn0_clear()

end subroutine writeemat
