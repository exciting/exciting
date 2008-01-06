
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
  use m_getunit
  use m_filedel
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'writeemat2'
  integer :: iq,un,qi

  ! initialise universal variables
  call init0
  call init1

  ! initialize q-point set
  call init2xs

  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)

!!$  ! k-point interval for process
!!$  kpari=firstofset(rank,nkpt)
!!$  kparf=lastofset(rank,nkpt)
  kpari=1
  kparf=nkpt

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

  ! resume task, second checkpoint index is q-point index
  qi=1
  call getunit(un)
  ! loop over q-points
  do iq = qi, nqpt
     ! call for q-point
     call ematq2(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): matrix elements of the &
          &exponentials finished for q-point:',iq
     call flushifc(unitout)
  end do

!!$  ! synchronize
!!$  call getunit(un)
!!$  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')
!!$
!!$  ! gather from processes
!!$  if ((procs.gt.1).and.(rank.eq.0)) call ematgather
!!$
!!$  ! synchronize
!!$  call getunit(un)
!!$  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

  write(unitout,'(a)') "Info("//trim(thisnam)//"): matrix elements of &
       &exponential expression finished"

  call findgntn0_clear

end subroutine writeemat2
