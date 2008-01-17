
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine df
  use modmain
  use modxs
  use modmpi
  use m_dfq
  use m_dfq2
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'df'
  integer :: iq

  call genfilname(setfilext=.true.)

  if (calledxs.eq.1) call init0

  ! initialise universal variables
  call init1

  ! save Gamma-point variables
  call tdsave0

  ! initialize q-point set
  call init2xs

  ! w-point parallelization for dielectric function
  partype='w'
  call genparidxran(partype)

  ! loop over q-points
  do iq=1,nqpt
     ! call for q-point
!!!     if (.not.gather) call dfq(iq)
     if (.not.gather) call dfq2(iq)
    write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn Sahm response &
          &function finished for q-point:',iq
  end do

  ! synchronize
  if (.not.gather) call barrier
  if ((procs.gt.1).and.(rank.eq.0)) call dfgather
  if (.not.gather) call barrier

  write(unitout,'(a)') "Info("//trim(thisnam)//"): Kohn-Sham response &
       &function finished"

  if (gather) then
     write(unitout,'(a)') "Info("//trim(thisnam)//"): gather option: &
          &exiting program"
     call xsfinit
     call terminate
  end if

end subroutine df
